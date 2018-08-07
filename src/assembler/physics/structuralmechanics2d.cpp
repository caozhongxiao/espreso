
#include "structuralmechanics2d.h"

#include "../step.h"
#include "../instance.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/evaluator/evaluator.h"
#include "../../basis/matrices/denseMatrix.h"

#include "../../config/ecf/physics/structuralmechanics.h"
#include "../../config/ecf/output.h"

#include "../../mesh/mesh.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/nodestore.h"
#include "../../mesh/store/fetidatastore.h"
#include "../../mesh/store/boundaryregionstore.h"
#include "../../mesh/store/elementsregionstore.h"
#include "../../mesh/store/surfacestore.h"

#include "../../solver/generic/SparseMatrix.h"
#include "../../solver/specific/sparsesolvers.h"

using namespace espreso;

StructuralMechanics2D::StructuralMechanics2D(Mesh *mesh, Instance *instance, Step *step, const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration)
: Physics("STRUCTURAL MECHANICS 2D", mesh, instance, step, &configuration, 2), StructuralMechanics(configuration, propertiesConfiguration, 2)
{

}

void StructuralMechanics2D::analyticRegularization(size_t domain, bool ortogonalCluster)
{
	if (_instance->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC";
	}

	Point center; size_t np; double norm;
	if (ortogonalCluster) {
		center = _cCenter[_mesh->elements->clusters[domain]];
		np = _cNp[_mesh->elements->clusters[domain]];
		norm = _cNorm[_mesh->elements->clusters[domain]].x;
	} else {
		center = _dCenter[domain];
		np = _dNp[domain];
		norm = _dNorm[domain].x;
	}

	_instance->N1[domain].rows = _instance->domainDOFCount[domain];
	_instance->N1[domain].cols = 3;
	_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
	_instance->N1[domain].type = 'G';

	_instance->N1[domain].dense_values.reserve(_instance->N1[domain].nnz);

	for (size_t c = 0; c < 2; c++) {
		std::vector<double> kernel = { 0, 0 };
		kernel[c] = 1 / std::sqrt(np);
		for (size_t i = 0; i < _instance->domainDOFCount[domain] / 2; i++) {
			_instance->N1[domain].dense_values.insert(_instance->N1[domain].dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < _mesh->nodes->dintervals[domain].size(); i++) {
		for (eslocal n = _mesh->nodes->dintervals[domain][i].begin; n < _mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = _mesh->nodes->coordinates->datatarray()[n] - center;
			_instance->N1[domain].dense_values.push_back(-p.y / norm);
			_instance->N1[domain].dense_values.push_back( p.x / norm);
		}
	}

	std::vector<eslocal> fixPoints;
	if (_BEMDomain[domain]) {
		fixPoints = std::vector<eslocal>(
				_mesh->FETIData->surfaceFixPoints.begin() + _mesh->FETIData->sFixPointsDistribution[domain],
				_mesh->FETIData->surfaceFixPoints.begin() + _mesh->FETIData->sFixPointsDistribution[domain + 1]);
	} else {
		fixPoints = std::vector<eslocal>(
				_mesh->FETIData->innerFixPoints.begin() + _mesh->FETIData->iFixPointsDistribution[domain],
				_mesh->FETIData->innerFixPoints.begin() + _mesh->FETIData->iFixPointsDistribution[domain + 1]);
	}

	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 3;
	Nt.cols = _instance->K[domain].cols;
	Nt.nnz  = 4 * fixPoints.size();
	Nt.type = 'G';

	std::vector<eslocal> &ROWS = Nt.CSR_I_row_indices;
	std::vector<eslocal> &COLS = Nt.CSR_J_col_indices;
	std::vector<double>  &VALS = Nt.CSR_V_values;

	ROWS.reserve(Nt.rows + 1);
	COLS.reserve(Nt.nnz);
	VALS.reserve(Nt.nnz);

	ROWS.push_back(1);
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());

	auto n2DOF = [&] (eslocal node) {
		auto dit = _mesh->nodes->dintervals[domain].begin();
		while (dit->end < node) { ++dit; }
		return dit->DOFOffset + node - dit->begin;
	};

	for (size_t c = 0; c < 2; c++) {
		for (size_t i = 0; i < fixPoints.size(); i++) {
			COLS.push_back(2 * n2DOF(fixPoints[i]) + c + 1);
		}
	}
	VALS.insert(VALS.end(), 2 * fixPoints.size(), 1);

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = _mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(2 * n2DOF(fixPoints[i]) + 0 + 1);
		COLS.push_back(2 * n2DOF(fixPoints[i]) + 1 + 1);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );

	_instance->RegMat[domain].MatMat(Nt, 'N', N);
	_instance->RegMat[domain].MatTranspose();
	_instance->RegMat[domain].RemoveLower();

	SparseSolverCPU NtN;
	NtN.ImportMatrix(_instance->RegMat[domain]);
	_instance->RegMat[domain].Clear();

	NtN.Factorization("Create RegMat");
	NtN.SolveMat_Sparse(Nt);
	NtN.Clear();

	_instance->RegMat[domain].MatMat(N, 'N', Nt);
	_instance->RegMat[domain].MatScale(_instance->K[domain].getDiagonalMaximum());
}

void StructuralMechanics2D::processBEM(eslocal domain, Matrices matrices)
{

}

void StructuralMechanics2D::assembleMaterialMatrix(eslocal node, const Point &p, const MaterialBaseConfiguration *mat, double temp, DenseMatrix &K) const
{
	double Ex, Ey, mi;

	switch (mat->linear_elastic_properties.model) {

	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
		Ex = Ey = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		mi = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		break;

	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
		ESINFO(ERROR) << "Implement ANISOTROPIC MATERIAL";
		break;

	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
		Ex = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		Ey = mat->linear_elastic_properties.young_modulus.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		mi = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		break;

	default:
		ESINFO(ERROR) << "Linear elasticity 2D not supports set material model";
	}

	switch (mat->linear_elastic_properties.model) {

	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
	{

		switch (_configuration.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
		{
			double k = Ex * (1 - mi) / ((1 + mi) * (1 - 2 * mi));
			K(node, 0) = k * 1;
			K(node, 1) = k * 1;
			K(node, 2) = k * ((1 - 2 * mi) / (2 * (1 - mi)));
			K(node, 3) = k * (mi / (1 - mi));
			K(node, 4) = 0;
			K(node, 5) = 0;
			K(node, 6) = k * (mi / (1 - mi));
			K(node, 7) = 0;
			K(node, 8) = 0;
			return;
		}

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
		{
			double k = Ex / (1 - mi * mi);
			K(node, 0) = k * 1;
			K(node, 1) = k * 1;
			K(node, 2) = k * ((1 -  mi) / 2);
			K(node, 3) = k * mi;
			K(node, 4) = 0;
			K(node, 5) = 0;
			K(node, 6) = k * mi;
			K(node, 7) = 0;
			K(node, 8) = 0;
			return;
		}

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
		{
			K.resize(K.rows(), 16);
			double k = Ex * (1 - mi) / ((1 + mi) * (1 - 2 * mi));
			K(node,  0) = k * 1;
			K(node,  1) = k * 1;
			K(node,  2) = k * ((1 - 2 * mi) / (2 * (1 - mi)));
			K(node,  3) = k * 1;

			K(node,  4) = k * (mi / (1 - mi));
			K(node,  5) = 0;
			K(node,  6) = k * (mi / (1 - mi));
			K(node,  7) = 0;
			K(node,  8) = k * (mi / (1 - mi));
			K(node,  9) = 0;

			K(node, 10) = k * (mi / (1 - mi));
			K(node, 11) = 0;
			K(node, 12) = 0;
			K(node, 13) = k * (mi / (1 - mi));
			K(node, 14) = k * (mi / (1 - mi));
			K(node, 15) = 0;
			return;
		}
		}
		break;
	}

	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
	{
		ESINFO(ERROR) << "IMPLEMENT: MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC";
		return;
	}

	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
	{
		ESINFO(ERROR) << "IMPLEMENT: MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC";
		return;
	}

	default:
		ESINFO(ERROR) << "Structural mechanics 2D not supports set material model";
	}
}

void StructuralMechanics2D::processElement(eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	auto nodes = _mesh->elements->nodes->cbegin() + eindex;
	auto epointer = _mesh->elements->epointers->datatarray()[eindex];
	const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[domain];
	const ECFExpressionVector *acceleration = NULL, *angular_velocity = NULL;
	Evaluator *thick = NULL;
	for (auto it = _configuration.load_steps_settings.at(_step->step + 1).acceleration.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).acceleration.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			acceleration = &it->second;
			break;
		}
	}
	for (auto it = _configuration.load_steps_settings.at(_step->step + 1).angular_velocity.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).angular_velocity.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			angular_velocity = &it->second;
			break;
		}
	}

	Evaluator *initial_temperature = NULL, *temperature = NULL;
	for (auto it = _configuration.load_steps_settings.at(_step->step + 1).temperature.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).temperature.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			temperature = it->second.evaluator;
			break;
		}
	}
	for (auto it = _configuration.initial_temperature.begin(); it != _configuration.initial_temperature.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			initial_temperature = it->second.evaluator;
			break;
		}
	}
	for (auto it = _configuration.thickness.begin(); it != _configuration.thickness.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			thick = it->second.evaluator;
			break;
		}
	}


	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);


	DenseMatrix Ce(4, 4), XY(1, 2), coordinates(nodes->size(), 2), J, invJ(2, 2), dND, B, precision, rhsT;
	DenseMatrix K(nodes->size(), 9), TE(nodes->size(), 2), thickness(nodes->size(), 1), inertia(nodes->size(), 2), dens(nodes->size(), 1), angvel(1, 3);
	DenseMatrix gpK(nodes->size(), 9), gpTE(1, 2), gpThickness(1, 1), gpInertia(1, 2), gpDens(1, 1);
	double detJ, temp, initTemp, CP = 1;
	Point center;

	const MaterialConfiguration* material = _mesh->materials[_mesh->elements->material->datatarray()[eindex]];

	for (size_t i = 0; i < nodes->size(); i++) {
		const Point &p = _mesh->nodes->coordinates->datatarray()[nodes->at(i)];
		if (initial_temperature != NULL) {
			initTemp = initial_temperature->evaluate(p, _step->currentTime, 0);
		}
		if (temperature != NULL) {
			temp = temperature->evaluate(p, _step->currentTime, 0);
		}
		if (acceleration != NULL) {
			inertia(i, 0) = acceleration->x.evaluator->evaluate(p, _step->currentTime, temp);
			inertia(i, 1) = acceleration->y.evaluator->evaluate(p, _step->currentTime, temp);
		}
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		dens(i, 0) = material->density.evaluator->evaluate(p, _step->currentTime, temp);

		thickness(i, 0) = thick != NULL ? thick->evaluate(p, _step->currentTime, temp) : 1;
		center += p;

		switch (material->linear_elastic_properties.model) {
		case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
			TE(i, 0) = TE(i, 1) = (temp - initTemp) * material->linear_elastic_properties.thermal_expansion.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
			break;
		case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
			TE(i, 0) = (temp - initTemp) * material->linear_elastic_properties.thermal_expansion.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
			TE(i, 1) = (temp - initTemp) * material->linear_elastic_properties.thermal_expansion.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Invalid LINEAR ELASTIC model.";
		}
		assembleMaterialMatrix(i, p, material, temp, K);
	}
	center /= nodes->size();
	if (angular_velocity != NULL) {
		angvel(0, 0) = acceleration->x.evaluator->evaluate(center, _step->currentTime, temp);
		angvel(0, 1) = acceleration->y.evaluator->evaluate(center, _step->currentTime, temp);
		angvel(0, 2) = acceleration->z.evaluator->evaluate(center, _step->currentTime, temp);
	}

	eslocal Ksize = 2 * nodes->size();

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if (matrices & (Matrices::K | Matrices::R)) {
		Ke.resize(Ksize, Ksize);
		Ke = 0;
	}
	if (matrices & Matrices::M) {
		Me.resize(Ksize, Ksize);
		Me = 0;
	}
	if (matrices & Matrices::R) {
		Re.resize(Ksize, 1);
		Re = 0;
	}
	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J.values());
		inverse2x2(J.values(), invJ.values(), detJ);

		gpThickness.multiply(N[gp], thickness);
		gpK.multiply(N[gp], K);
		dND.multiply(invJ, dN[gp]);
		gpDens.multiply(N[gp], dens);

		if (matrices & Matrices::f) {
			gpTE.multiply(N[gp], TE);
			gpInertia.multiply(N[gp], inertia);
			XY.multiply(N[gp], coordinates);
		}

		if (matrices & Matrices::M) {
			Me.multiply(N[gp], N[gp], gpDens(0, 0) * detJ * weighFactor[gp] * CP, 1, true);
		}

		switch (_configuration.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:

			Ce.resize(3, 3);
			Ce(0, 0) = gpK(0, 0);
			Ce(1, 1) = gpK(0, 1);
			Ce(2, 2) = gpK(0, 2);
			Ce(0, 1) = gpK(0, 3);
			Ce(0, 2) = gpK(0, 4);
			Ce(1, 2) = gpK(0, 5);
			Ce(1, 0) = gpK(0, 6);
			Ce(2, 0) = gpK(0, 7);
			Ce(2, 1) = gpK(0, 8);

			B.resize(Ce.rows(), Ksize);
			distribute3x2(B.values(), dND.values(), dND.rows(), dND.columns());

			if (matrices & (Matrices::K | Matrices::R)) {
				Ke.multiply(B, Ce * B, detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
			}

			if (matrices & Matrices::f) {
				precision.resize(Ce.rows(), 1);
				precision(0, 0) = gpTE(0, 1);
				precision(1, 0) = gpTE(0, 1);
				precision(2, 0) = 0;
				rhsT.multiply(B, Ce * precision, detJ * weighFactor[gp] * gpThickness(0, 0), 0, true, false);
				for (eslocal i = 0; i < Ksize; i++) {
					fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % nodes->size()) * gpInertia(0, i / nodes->size());
					fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % nodes->size()) * XY(0, i / nodes->size()) * pow(angvel(0, 2), 2);
					fe(i, 0) += rhsT(i, 0);
				}
			}
			break;

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:

			Ce.resize(4, 4);
			Ce(0, 0) = gpK(0,  0);
			Ce(1, 1) = gpK(0,  1);
			Ce(2, 2) = gpK(0,  2);
			Ce(3, 3) = gpK(0,  3);
			Ce(0, 1) = gpK(0,  4);
			Ce(0, 2) = gpK(0,  5);
			Ce(0, 3) = gpK(0,  6);
			Ce(1, 2) = gpK(0,  7);
			Ce(1, 3) = gpK(0,  8);
			Ce(2, 3) = gpK(0,  9);
			Ce(1, 0) = gpK(0, 10);
			Ce(2, 0) = gpK(0, 11);
			Ce(2, 1) = gpK(0, 12);
			Ce(3, 0) = gpK(0, 13);
			Ce(3, 1) = gpK(0, 14);
			Ce(3, 2) = gpK(0, 15);

			B.resize(Ce.rows(), Ksize);
			distribute4x2(B.values(), dND.values(), dND.rows(), dND.columns());
			for(size_t i = 0; i < N[gp].columns(); i++) {
				B(2, i) = N[gp](0, i) / XY(0, 0);
			}

			if (matrices & (Matrices::K | Matrices::R)) {
				Ke.multiply(B, Ce * B, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 1, true);
			}

			if (matrices & Matrices::f) {
				precision.resize(Ce.rows(), 1);
				precision(0, 0) = gpTE(0, 0);
				precision(1, 0) = gpTE(0, 1);
				precision(2, 0) = precision(3, 0) = 0;
				rhsT.multiply(B, Ce * precision, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 0, true, false);
				for (eslocal i = 0; i < Ksize; i++) {
					fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % nodes->size()) * gpInertia(0, i / nodes->size());
					fe(i, 0) += rhsT(i, 0);
				}
				for (eslocal i = 0; i < Ksize / 2; i++) {
					fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % nodes->size()) * XY(0, 0) * pow(angvel(0, 1), 2);
					fe(Ksize / 2 + i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % nodes->size()) * XY(0, 1) * pow(angvel(0, 1), 2);
				}
			}
			break;
		}
	}
}

void StructuralMechanics2D::processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	ESINFO(ERROR) << "Structural mechanics 2D cannot process face";
}

void StructuralMechanics2D::processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	const Evaluator *pressure = NULL, *thickness = NULL;
	auto it = _configuration.load_steps_settings.at(_step->step + 1).normal_pressure.find(region->name);
	if (it != _configuration.load_steps_settings.at(_step->step + 1).normal_pressure.end()) {
		pressure = it->second.evaluator;
	}
	for (auto it = _configuration.thickness.begin(); it != _configuration.thickness.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			thickness = it->second.evaluator;
			break;
		}
	}

	if (pressure == NULL) {
		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}
	if (!(matrices & (Matrices::K | Matrices::f))) {
		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}

	auto nodes = region->elements->cbegin() + eindex;
	auto epointer = region->epointers->datatarray()[eindex];
	const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[domain];

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	DenseMatrix coordinates(nodes->size(), 2), dND(1, 2), P(nodes->size(), 1), normal(2, 1), matThickness(nodes->size(), 1), XY(1, 2);
	DenseMatrix gpP(1, 1), gpQ(1, 2), gpThickness(1, 1);

	eslocal Ksize = 2 * nodes->size();
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);

	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	for (size_t n = 0; n < nodes->size(); n++) {
		double temp = 0;
		const Point &p = _mesh->nodes->coordinates->datatarray()[nodes->at(n)];
		coordinates(n, 0) = p.x;
		coordinates(n, 1) = p.y;
		if (pressure != NULL) {
			P(n, 0) = pressure->evaluate(p, _step->currentTime, temp);
		}
		matThickness(n, 0) = thickness != NULL ? thickness->evaluate(p, _step->currentTime, temp) : 1;
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		Point n(-dND(0, 1), dND(0, 0), 0);
		// e->rotateOutside(e->parentElements()[0], _mesh->coordinates(), n);
		normal(0, 0) = n.x / J;
		normal(1, 0) = n.y / J;
		gpP.multiply(N[gp], P);
		gpQ.multiply(normal, gpP);
		gpThickness.multiply(N[gp], matThickness);

		switch (_configuration.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			for (eslocal i = 0; i < Ksize; i++) {
				fe(i, 0) += gpThickness(0, 0) * J * weighFactor[gp] * N[gp](0, i % nodes->size()) * gpQ(0, i / nodes->size());
			}
			break;

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			XY.multiply(N[gp], coordinates);
			for (eslocal i = 0; i < Ksize; i++) {
				fe(i, 0) += gpThickness(0, 0) * J * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % nodes->size()) * gpQ(0, i / nodes->size());
			}
			break;
		}
	}
}

void StructuralMechanics2D::processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
//	if (
//			e->hasProperty(Property::FORCE_X, step.step) ||
//			e->hasProperty(Property::FORCE_Y, step.step)) {
//
//		Ke.resize(0, 0);
//		Me.resize(0, 0);
//		Re.resize(0, 0);
//		fe.resize(pointDOFs().size(), 0);
//
////		fe(0, 0) = e->sumProperty(Property::FORCE_X, step.step, _mesh->coordinates()[e->node(0)], step.currentTime, 0, 0);
////		fe(1, 0) = e->sumProperty(Property::FORCE_Y, step.step, _mesh->coordinates()[e->node(0)], step.currentTime, 0, 0);
//		return;
//	}
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
}

void StructuralMechanics2D::postProcessElement(eslocal eindex)
{

}

void StructuralMechanics2D::processSolution()
{
}





