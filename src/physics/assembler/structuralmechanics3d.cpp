
#include "../assembler/structuralmechanics3d.h"

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

StructuralMechanics3D::StructuralMechanics3D(Mesh *mesh, Instance *instance, Step *step, const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration)
: Physics("STRUCTURAL MECHANICS 3D", mesh, instance, step, &configuration, 3), StructuralMechanics(configuration, propertiesConfiguration, 3)
{

}

void StructuralMechanics3D::analyticRegularization(size_t domain, bool ortogonalCluster)
{
	if (_instance->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC";
	}

	Point center = _dCenter[domain], norm = _dNorm[domain];
	double r44 = _dr44[domain], r45 = _dr45[domain], r46 = _dr46[domain], r55 = _dr55[domain], r56 = _dr56[domain];
	size_t np = _dNp[domain];

	if (ortogonalCluster) {
		size_t cluster = _mesh->elements->clusters[domain];
		center = _cCenter[cluster], norm = _cNorm[cluster];
		r44 = _cr44[cluster], r45 = _cr45[cluster], r46 = _cr46[cluster], r55 = _cr55[cluster], r56 = _cr56[cluster];
		np = _cNp[cluster];
	} else {
		center = _dCenter[domain], norm = _dNorm[domain];
		r44 = _dr44[domain], r45 = _dr45[domain], r46 = _dr46[domain], r55 = _dr55[domain], r56 = _dr56[domain];
		np = _dNp[domain];
	}

	_instance->N1[domain].rows = _instance->domainDOFCount[domain];
	_instance->N1[domain].cols = 6;
	_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
	_instance->N1[domain].type = 'G';

	_instance->N1[domain].dense_values.reserve(_instance->N1[domain].nnz);

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };
		kernel[c] = 1 / std::sqrt(np);
		for (size_t i = 0; i < _instance->domainDOFCount[domain] / 3; i++) {
			_instance->N1[domain].dense_values.insert(_instance->N1[domain].dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < _mesh->nodes->dintervals[domain].size(); i++) {
		for (eslocal n = _mesh->nodes->dintervals[domain][i].begin; n < _mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = _mesh->nodes->coordinates->datatarray()[n] - center;
			_instance->N1[domain].dense_values.push_back(-p.y / norm.x);
			_instance->N1[domain].dense_values.push_back( p.x / norm.x);
			_instance->N1[domain].dense_values.push_back(             0);
		}
	}

	for (size_t i = 0; i < _mesh->nodes->dintervals[domain].size(); i++) {
		for (eslocal n = _mesh->nodes->dintervals[domain][i].begin; n < _mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = _mesh->nodes->coordinates->datatarray()[n] - center;
			_instance->N1[domain].dense_values.push_back((-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y);
			_instance->N1[domain].dense_values.push_back((   0 - r45 / r44 * ( p.x / norm.x)) / norm.y);
			_instance->N1[domain].dense_values.push_back(( p.x - r45 / r44 * (   0 / norm.x)) / norm.y);
		}
	}

	for (size_t i = 0; i < _mesh->nodes->dintervals[domain].size(); i++) {
		for (eslocal n = _mesh->nodes->dintervals[domain][i].begin; n < _mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = _mesh->nodes->coordinates->datatarray()[n] - center;
			_instance->N1[domain].dense_values.push_back((   0 - r56 / r55 * ((-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y) - r46 / r44 * (-p.y / norm.x)) / norm.z);
			_instance->N1[domain].dense_values.push_back((-p.z - r56 / r55 * ((   0 - r45 / r44 * ( p.x / norm.x)) / norm.y) - r46 / r44 * ( p.x / norm.x)) / norm.z);
			_instance->N1[domain].dense_values.push_back(( p.y - r56 / r55 * (( p.x - r45 / r44 * (   0 / norm.x)) / norm.y) - r46 / r44 * (   0 / norm.x)) / norm.z);
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
	Nt.rows = 6;
	Nt.cols = _instance->K[domain].cols;
	Nt.nnz  = 9 * fixPoints.size();
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
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());

	auto n2DOF = [&] (eslocal node) {
		auto dit = _mesh->nodes->dintervals[domain].begin();
		while (dit->end < node) { ++dit; }
		return dit->DOFOffset + node - dit->begin;
	};

	for (size_t c = 0; c < 3; c++) {
		for (size_t i = 0; i < fixPoints.size(); i++) {
			COLS.push_back(3 * n2DOF(fixPoints[i]) + c + 1);
		}
	}
	VALS.insert(VALS.end(), 3 * fixPoints.size(), 1);

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = _mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 0 + 1);
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 1 + 1);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = _mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 0 + 1);
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 2 + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = _mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 1 + 1);
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 2 + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.y);
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

void StructuralMechanics3D::processBEM(eslocal domain, Matrices matrices)
{

}

void StructuralMechanics3D::assembleMaterialMatrix(eslocal node, const Point &p, const MaterialBaseConfiguration *mat, double temp, DenseMatrix &K) const
{
	double Ex, Ey, Ez, miXY, miXZ, miYZ, Gx, Gy, Gz;

	switch (mat->linear_elastic_properties.model) {

	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC: {
		Ex = Ey = Ez = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		miXY = miXZ = miYZ = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);

		double EE = Ex / ((1 + miXY) * (1 - 2 * miXY));

		K(node,  0) = EE * (1.0 - miXY);
		K(node,  1) = EE * (1.0 - miXY);
		K(node,  2) = EE * (1.0 - miXY);
		K(node,  3) = EE * (0.5 - miXY);
		K(node,  4) = EE * (0.5 - miXY);
		K(node,  5) = EE * (0.5 - miXY);

		K(node,  6) = EE * miXY;
		K(node,  7) = EE * miXY;
		K(node,  8) = 0;
		K(node,  9) = 0;
		K(node, 10) = 0;
		K(node, 11) = EE * miXY;
		K(node, 12) = 0;
		K(node, 13) = 0;
		K(node, 14) = 0;
		K(node, 15) = 0;
		K(node, 16) = 0;
		K(node, 17) = 0;
		K(node, 18) = 0;
		K(node, 19) = 0;
		K(node, 20) = 0;
		K(node, 21) = EE * miXY;
		K(node, 22) = EE * miXY;
		K(node, 23) = EE * miXY;
		K(node, 24) = 0;
		K(node, 25) = 0;
		K(node, 26) = 0;
		K(node, 27) = 0;
		K(node, 28) = 0;
		K(node, 29) = 0;
		K(node, 30) = 0;
		K(node, 31) = 0;
		K(node, 32) = 0;
		K(node, 33) = 0;
		K(node, 34) = 0;
		K(node, 35) = 0;
	} break;

	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC: {
		K(node,  0) = mat->linear_elastic_properties.anisotropic.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		K(node,  1) = mat->linear_elastic_properties.anisotropic.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		K(node,  2) = mat->linear_elastic_properties.anisotropic.get(2, 2).evaluator->evaluate(p, _step->currentTime, temp);
		K(node,  3) = mat->linear_elastic_properties.anisotropic.get(3, 3).evaluator->evaluate(p, _step->currentTime, temp);
		K(node,  4) = mat->linear_elastic_properties.anisotropic.get(4, 4).evaluator->evaluate(p, _step->currentTime, temp);
		K(node,  5) = mat->linear_elastic_properties.anisotropic.get(5, 5).evaluator->evaluate(p, _step->currentTime, temp);

		K(node,  6) = mat->linear_elastic_properties.anisotropic.get(0, 1).evaluator->evaluate(p, _step->currentTime, temp);
		K(node,  7) = mat->linear_elastic_properties.anisotropic.get(0, 2).evaluator->evaluate(p, _step->currentTime, temp);
		K(node,  8) = mat->linear_elastic_properties.anisotropic.get(0, 3).evaluator->evaluate(p, _step->currentTime, temp);
		K(node,  9) = mat->linear_elastic_properties.anisotropic.get(0, 4).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 10) = mat->linear_elastic_properties.anisotropic.get(0, 5).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 11) = mat->linear_elastic_properties.anisotropic.get(1, 2).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 12) = mat->linear_elastic_properties.anisotropic.get(1, 3).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 13) = mat->linear_elastic_properties.anisotropic.get(1, 4).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 14) = mat->linear_elastic_properties.anisotropic.get(1, 5).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 15) = mat->linear_elastic_properties.anisotropic.get(2, 3).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 16) = mat->linear_elastic_properties.anisotropic.get(2, 4).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 17) = mat->linear_elastic_properties.anisotropic.get(2, 5).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 18) = mat->linear_elastic_properties.anisotropic.get(3, 4).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 19) = mat->linear_elastic_properties.anisotropic.get(3, 5).evaluator->evaluate(p, _step->currentTime, temp);
		K(node, 20) = mat->linear_elastic_properties.anisotropic.get(4, 5).evaluator->evaluate(p, _step->currentTime, temp);

		K(node, 21) = K(node,  6);
		K(node, 22) = K(node,  7);
		K(node, 23) = K(node, 11);
		K(node, 24) = K(node,  8);
		K(node, 25) = K(node, 12);
		K(node, 26) = K(node, 15);
		K(node, 27) = K(node,  9);
		K(node, 28) = K(node, 13);
		K(node, 29) = K(node, 16);
		K(node, 30) = K(node, 18);
		K(node, 31) = K(node, 10);
		K(node, 32) = K(node, 14);
		K(node, 33) = K(node, 17);
		K(node, 34) = K(node, 19);
		K(node, 35) = K(node, 20);
	} break;

	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC: {
		Ex = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		Ey = mat->linear_elastic_properties.young_modulus.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		Ez = mat->linear_elastic_properties.young_modulus.get(2, 2).evaluator->evaluate(p, _step->currentTime, temp);

		miXY = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		miXZ = mat->linear_elastic_properties.poisson_ratio.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		miYZ = mat->linear_elastic_properties.poisson_ratio.get(2, 2).evaluator->evaluate(p, _step->currentTime, temp);

		Gx = mat->linear_elastic_properties.shear_modulus.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		Gy = mat->linear_elastic_properties.shear_modulus.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		Gz = mat->linear_elastic_properties.shear_modulus.get(2, 2).evaluator->evaluate(p, _step->currentTime, temp);

		double miYX = miXY * Ey / Ex;
		double miZY = miYZ * Ez / Ey;
		double miZX = miXZ * Ex / Ez;

		double ksi = 1 - (miXY * miYX + miYZ * miZY + miXZ * miZX) - (miXY * miYZ * miZX + miYX * miZY * miXZ);

		double dxx = Ex * (1 - miYZ * miZY) / ksi;
		double dxy = Ey * (miXY + miXZ * miZY) / ksi;
		double dxz = Ez * (miXZ + miYZ * miXY)  /ksi;
		double dyy = Ey * (1 - miXZ * miZX) / ksi;
		double dyz = Ez * (miYZ + miYX * miXZ) / ksi;
		double dzz = Ez * (1 - miYX * miXY) / ksi;

		K(node,  0) = dxx;
		K(node,  1) = dyy;
		K(node,  2) = dzz;
		K(node,  3) = Gx;
		K(node,  4) = Gz;
		K(node,  5) = Gy;

		K(node,  6) = dxy;
		K(node,  7) = dxz;
		K(node,  8) = 0;
		K(node,  9) = 0;
		K(node, 10) = 0;
		K(node, 11) = dyz;
		K(node, 12) = 0;
		K(node, 13) = 0;
		K(node, 14) = 0;
		K(node, 15) = 0;
		K(node, 16) = 0;
		K(node, 17) = 0;
		K(node, 18) = 0;
		K(node, 19) = 0;
		K(node, 20) = 0;
		K(node, 21) = dxy;
		K(node, 22) = dxz;
		K(node, 23) = dyz;
		K(node, 24) = 0;
		K(node, 25) = 0;
		K(node, 26) = 0;
		K(node, 27) = 0;
		K(node, 28) = 0;
		K(node, 29) = 0;
		K(node, 30) = 0;
		K(node, 31) = 0;
		K(node, 32) = 0;
		K(node, 33) = 0;
		K(node, 34) = 0;
		K(node, 35) = 0;
	} break;

	default:
		ESINFO(ERROR) << "Structural mechanics 3D not supports set material model";
	}
}

void StructuralMechanics3D::processElement(eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	auto nodes = _mesh->elements->procNodes->cbegin() + eindex;
	auto epointer = _mesh->elements->epointers->datatarray()[eindex];
	const ECFExpressionVector *acceleration = NULL;
	for (auto it = _configuration.load_steps_settings.at(_step->step + 1).acceleration.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).acceleration.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			acceleration = &it->second;
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


	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	DenseMatrix Ce(6, 6), coordinates(nodes->size(), 3), J, invJ(3, 3), dND, B, precision, rhsT;
	DenseMatrix K(nodes->size(), 36), TE(nodes->size(), 3), inertia(nodes->size(), 3), dens(nodes->size(), 1);
	DenseMatrix gpK(nodes->size(), 36), gpTE(1, 3), gpInertia(1, 3), gpDens(1, 1);
	double detJ, temp = 275.15, initTemp = 275.15, CP = 1;

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
			inertia(i, 2) = acceleration->z.evaluator->evaluate(p, _step->currentTime, temp);
		}
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		coordinates(i, 2) = p.z;
		dens(i, 0) = material->density.evaluator->evaluate(p, _step->currentTime, temp);
		switch (material->linear_elastic_properties.model) {
		case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
			TE(i, 0) = TE(i, 1) = TE(i, 2) = (temp - initTemp) * material->linear_elastic_properties.thermal_expansion.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
			break;
		case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
		case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
			TE(i, 0) = (temp - initTemp) * material->linear_elastic_properties.thermal_expansion.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
			TE(i, 1) = (temp - initTemp) * material->linear_elastic_properties.thermal_expansion.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
			TE(i, 2) = (temp - initTemp) * material->linear_elastic_properties.thermal_expansion.get(2, 2).evaluator->evaluate(p, _step->currentTime, temp);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Invalid LINEAR ELASTIC model.";
		}
		assembleMaterialMatrix(i, p, material, temp, K);
	}

	eslocal Ksize = 3 * nodes->size();

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
		detJ = determinant3x3(J.values());
		if (detJ <= 0) {
			printInvalidElement(eindex);
			ESINFO(ERROR) << "Invalid element detected - check input data.";
		}
		inverse3x3(J.values(), invJ.values(), detJ);

		gpK.multiply(N[gp], K);
		dND.multiply(invJ, dN[gp]);
		gpDens.multiply(N[gp], dens);

		if (matrices & Matrices::f) {
			gpTE.multiply(N[gp], TE);
			gpInertia.multiply(N[gp], inertia);
		}

		if (matrices & Matrices::M) {
			Me.multiply(N[gp], N[gp], gpDens(0, 0) * detJ * weighFactor[gp] * CP, 1, true);
		}

		Ce.resize(6, 6);
		size_t k = 0;
		for (size_t i = 0; i < 6; i++) {
			Ce(i, i) = gpK(0, k++);
		}
		for (size_t i = 0; i < 6; i++) {
			for (size_t j = i + 1; j < 6; j++) {
				Ce(i, j) = gpK(0, k++);
			}
		}
		for (size_t i = 0; i < 6; i++) {
			for (size_t j = 0; j < i; j++) {
				Ce(i, j) = gpK(0, k++);
			}
		}

		B.resize(Ce.rows(), Ksize);
		distribute6x3(B.values(), dND.values(), dND.rows(), dND.columns());

		if (matrices & Matrices::K) {
			Ke.multiply(B, Ce * B, detJ * weighFactor[gp], 1, true);
		}

		if (matrices & Matrices::f) {
			precision.resize(Ce.rows(), 1);
			precision(0, 0) = gpTE(0, 0);
			precision(1, 0) = gpTE(0, 1);
			precision(2, 0) = gpTE(0, 2);
			precision(3, 0) = precision(4, 0) = precision(5, 0) = 0;

			rhsT.multiply(B, Ce * precision, detJ * weighFactor[gp], 0, true, false);
			for (eslocal i = 0; i < Ksize; i++) {
				fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * N[gp](0, i % nodes->size()) * gpInertia(0, i / nodes->size());
				fe(i, 0) += rhsT(i, 0);
			}
		}
	}
}

void StructuralMechanics3D::processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	const Evaluator *pressure = NULL;
	auto it = _configuration.load_steps_settings.at(_step->step + 1).normal_pressure.find(region->name);
	if (it != _configuration.load_steps_settings.at(_step->step + 1).normal_pressure.end()) {
		pressure = it->second.evaluator;
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

	auto nodes = region->procNodes->cbegin() + findex;
	auto epointer = region->epointers->datatarray()[findex];

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	DenseMatrix coordinates(nodes->size(), 3), dND(1, 3), P(nodes->size(), 1), normal(1, 3);
	DenseMatrix gpP(1, 1), gpQ(1, 3);

	eslocal Ksize = 3 * nodes->size();
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);

	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	for (size_t n = 0; n < nodes->size(); n++) {
		const Point &p = _mesh->nodes->coordinates->datatarray()[nodes->at(n)];
		coordinates(n, 0) = p.x;
		coordinates(n, 1) = p.y;
		coordinates(n, 2) = p.z;
		P(n, 0) += pressure->evaluate(p, 0, _step->currentTime);
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		// e->rotateOutside(e->parentElements()[0], _mesh->coordinates(), va);
		double J = va.norm();
		normal(0, 0) = va.x / va.norm();
		normal(0, 1) = va.y / va.norm();
		normal(0, 2) = va.z / va.norm();

		gpP.multiply(N[gp], P);
		gpQ.multiply(normal, gpP, 1, 0, true);

		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % nodes->size()) * gpQ(0, i / nodes->size());
		}
	}
}

void StructuralMechanics3D::processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
}

void StructuralMechanics3D::processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
//	if (
//			e->hasProperty(Property::FORCE_X, step.step) ||
//			e->hasProperty(Property::FORCE_Y, step.step) ||
//			e->hasProperty(Property::FORCE_Z, step.step)) {
//
//		Ke.resize(0, 0);
//		Me.resize(0, 0);
//		Re.resize(0, 0);
//		fe.resize(pointDOFs().size(), 0);
//
////		fe(0, 0) = e->sumProperty(Property::FORCE_X, step.step, _mesh->coordinates()[e->node(0)], step.currentTime, 0, 0);
////		fe(1, 0) = e->sumProperty(Property::FORCE_Y, step.step, _mesh->coordinates()[e->node(0)], step.currentTime, 0, 0);
////		fe(2, 0) = e->sumProperty(Property::FORCE_Z, step.step, _mesh->coordinates()[e->node(0)], step.currentTime, 0, 0);
//		return;
//	}
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
}

void StructuralMechanics3D::postProcessElement(eslocal eindex)
{

}

void StructuralMechanics3D::processSolution()
{
}





