
#include "structuralmechanics3d.controller.h"
#include "../kernels/structuralmechanics3d.kernel.h"

#include "../../step.h"
#include "../../instance.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../config/ecf/environment.h"
#include "../../../config/ecf/physics/structuralmechanics.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/nodestore.h"

using namespace espreso;

StructuralMechanics3DControler::StructuralMechanics3DControler(
		Mesh &mesh, const Step &step,
		const StructuralMechanicsGlobalSettings &gSettings,
		const StructuralMechanicsStepSettings &sSettings,
		const StructuralMechanicsOutputSettings &oSettings)
: StructuralMechanicsControler(mesh, step, gSettings, sSettings, oSettings)
{
	_kernel = new StructuralMechanics3DKernel(_globalSettings, _outputSettings);

	Point defaultMotion(0, 0, 0);
	double defaultHeat = 0; //273.15;

	_ntemperature.data = new serializededata<eslocal, double>(1, _nDistribution);
	_ncoordinate.data = new serializededata<eslocal, double>(3, _nDistribution);

	_boundaries.resize(_mesh.boundaryRegions.size());
}

StructuralMechanics3DControler::~StructuralMechanics3DControler()
{
	delete _kernel;
}

void StructuralMechanics3DControler::analyticRegularization(size_t domain, bool ortogonalCluster)
{
//	if (_instance->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
//		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC";
//	}
//
//	Point center = _dCenter[domain], norm = _dNorm[domain];
//	double r44 = _dr44[domain], r45 = _dr45[domain], r46 = _dr46[domain], r55 = _dr55[domain], r56 = _dr56[domain];
//	size_t np = _dNp[domain];
//
//	if (ortogonalCluster) {
//		size_t cluster = _mesh->elements->clusters[domain];
//		center = _cCenter[cluster], norm = _cNorm[cluster];
//		r44 = _cr44[cluster], r45 = _cr45[cluster], r46 = _cr46[cluster], r55 = _cr55[cluster], r56 = _cr56[cluster];
//		np = _cNp[cluster];
//	} else {
//		center = _dCenter[domain], norm = _dNorm[domain];
//		r44 = _dr44[domain], r45 = _dr45[domain], r46 = _dr46[domain], r55 = _dr55[domain], r56 = _dr56[domain];
//		np = _dNp[domain];
//	}
//
//	_instance->N1[domain].rows = _instance->domainDOFCount[domain];
//	_instance->N1[domain].cols = 6;
//	_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
//	_instance->N1[domain].type = 'G';
//
//	_instance->N1[domain].dense_values.reserve(_instance->N1[domain].nnz);
//
//	for (size_t c = 0; c < 3; c++) {
//		std::vector<double> kernel = { 0, 0, 0 };
//		kernel[c] = 1 / std::sqrt(np);
//		for (size_t i = 0; i < _instance->domainDOFCount[domain] / 3; i++) {
//			_instance->N1[domain].dense_values.insert(_instance->N1[domain].dense_values.end(), kernel.begin(), kernel.end());
//		}
//	}
//
//	for (size_t i = 0; i < _mesh->nodes->dintervals[domain].size(); i++) {
//		for (eslocal n = _mesh->nodes->dintervals[domain][i].begin; n < _mesh->nodes->dintervals[domain][i].end; ++n) {
//			Point p = _mesh->nodes->coordinates->datatarray()[n] - center;
//			_instance->N1[domain].dense_values.push_back(-p.y / norm.x);
//			_instance->N1[domain].dense_values.push_back( p.x / norm.x);
//			_instance->N1[domain].dense_values.push_back(             0);
//		}
//	}
//
//	for (size_t i = 0; i < _mesh->nodes->dintervals[domain].size(); i++) {
//		for (eslocal n = _mesh->nodes->dintervals[domain][i].begin; n < _mesh->nodes->dintervals[domain][i].end; ++n) {
//			Point p = _mesh->nodes->coordinates->datatarray()[n] - center;
//			_instance->N1[domain].dense_values.push_back((-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y);
//			_instance->N1[domain].dense_values.push_back((   0 - r45 / r44 * ( p.x / norm.x)) / norm.y);
//			_instance->N1[domain].dense_values.push_back(( p.x - r45 / r44 * (   0 / norm.x)) / norm.y);
//		}
//	}
//
//	for (size_t i = 0; i < _mesh->nodes->dintervals[domain].size(); i++) {
//		for (eslocal n = _mesh->nodes->dintervals[domain][i].begin; n < _mesh->nodes->dintervals[domain][i].end; ++n) {
//			Point p = _mesh->nodes->coordinates->datatarray()[n] - center;
//			_instance->N1[domain].dense_values.push_back((   0 - r56 / r55 * ((-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y) - r46 / r44 * (-p.y / norm.x)) / norm.z);
//			_instance->N1[domain].dense_values.push_back((-p.z - r56 / r55 * ((   0 - r45 / r44 * ( p.x / norm.x)) / norm.y) - r46 / r44 * ( p.x / norm.x)) / norm.z);
//			_instance->N1[domain].dense_values.push_back(( p.y - r56 / r55 * (( p.x - r45 / r44 * (   0 / norm.x)) / norm.y) - r46 / r44 * (   0 / norm.x)) / norm.z);
//		}
//	}
//
//	std::vector<eslocal> fixPoints;
//	if (_BEMDomain[domain]) {
//		fixPoints = std::vector<eslocal>(
//				_mesh->FETIData->surfaceFixPoints.begin() + _mesh->FETIData->sFixPointsDistribution[domain],
//				_mesh->FETIData->surfaceFixPoints.begin() + _mesh->FETIData->sFixPointsDistribution[domain + 1]);
//	} else {
//		fixPoints = std::vector<eslocal>(
//				_mesh->FETIData->innerFixPoints.begin() + _mesh->FETIData->iFixPointsDistribution[domain],
//				_mesh->FETIData->innerFixPoints.begin() + _mesh->FETIData->iFixPointsDistribution[domain + 1]);
//	}
//
//	SparseMatrix Nt; // CSR matice s DOFY
//	Nt.rows = 6;
//	Nt.cols = _instance->K[domain].cols;
//	Nt.nnz  = 9 * fixPoints.size();
//	Nt.type = 'G';
//
//	std::vector<eslocal> &ROWS = Nt.CSR_I_row_indices;
//	std::vector<eslocal> &COLS = Nt.CSR_J_col_indices;
//	std::vector<double>  &VALS = Nt.CSR_V_values;
//
//	ROWS.reserve(Nt.rows + 1);
//	COLS.reserve(Nt.nnz);
//	VALS.reserve(Nt.nnz);
//
//	ROWS.push_back(1);
//	ROWS.push_back(ROWS.back() + fixPoints.size());
//	ROWS.push_back(ROWS.back() + fixPoints.size());
//	ROWS.push_back(ROWS.back() + fixPoints.size());
//	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
//	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
//	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
//
//	auto n2DOF = [&] (eslocal node) {
//		auto dit = _mesh->nodes->dintervals[domain].begin();
//		while (dit->end < node) { ++dit; }
//		return dit->DOFOffset + node - dit->begin;
//	};
//
//	for (size_t c = 0; c < 3; c++) {
//		for (size_t i = 0; i < fixPoints.size(); i++) {
//			COLS.push_back(3 * n2DOF(fixPoints[i]) + c + 1);
//		}
//	}
//	VALS.insert(VALS.end(), 3 * fixPoints.size(), 1);
//
//	for (size_t i = 0; i < fixPoints.size(); i++) {
//		const Point &p = _mesh->nodes->coordinates->datatarray()[fixPoints[i]];
//		COLS.push_back(3 * n2DOF(fixPoints[i]) + 0 + 1);
//		COLS.push_back(3 * n2DOF(fixPoints[i]) + 1 + 1);
//		VALS.push_back(-p.y);
//		VALS.push_back( p.x);
//	}
//
//	for (size_t i = 0; i < fixPoints.size(); i++) {
//		const Point &p = _mesh->nodes->coordinates->datatarray()[fixPoints[i]];
//		COLS.push_back(3 * n2DOF(fixPoints[i]) + 0 + 1);
//		COLS.push_back(3 * n2DOF(fixPoints[i]) + 2 + 1);
//		VALS.push_back(-p.z);
//		VALS.push_back( p.x);
//	}
//
//	for (size_t i = 0; i < fixPoints.size(); i++) {
//		const Point &p = _mesh->nodes->coordinates->datatarray()[fixPoints[i]];
//		COLS.push_back(3 * n2DOF(fixPoints[i]) + 1 + 1);
//		COLS.push_back(3 * n2DOF(fixPoints[i]) + 2 + 1);
//		VALS.push_back(-p.z);
//		VALS.push_back( p.y);
//	}
//
//	SparseMatrix N;
//	Nt.MatTranspose( N );
//
//	_instance->RegMat[domain].MatMat(Nt, 'N', N);
//	_instance->RegMat[domain].MatTranspose();
//	_instance->RegMat[domain].RemoveLower();
//
//	SparseSolverCPU NtN;
//	NtN.ImportMatrix(_instance->RegMat[domain]);
//	_instance->RegMat[domain].Clear();
//
//	NtN.Factorization("Create RegMat");
//	NtN.SolveMat_Sparse(Nt);
//	NtN.Clear();
//
//	_instance->RegMat[domain].MatMat(N, 'N', Nt);
//	_instance->RegMat[domain].MatScale(_instance->K[domain].getDiagonalMaximum());
}

void StructuralMechanics3DControler::dirichletIndices(std::vector<std::vector<eslocal> > &indices)
{
	indices.resize(3);

	for (auto it = _stepSettings.displacement.regions.begin(); it != _stepSettings.displacement.regions.end(); ++it) {
		BoundaryRegionStore *region = _mesh.bregion(it->first);
		indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
	}

	for (auto it = _stepSettings.displacement.intersections.begin(); it != _stepSettings.displacement.intersections.end(); ++it) {
		BoundaryRegionsIntersectionStore *region = _mesh.ibregion(it->first);
		indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
	}
	_dirichletSize = indices[0].size();
}

void StructuralMechanics3DControler::dirichletValues(std::vector<double> &values)
{
	values.resize(_dirichletSize);

//	size_t offset = 0;
//	for (auto it = _stepSettings.displacement.regions.begin(); it != _stepSettings.displacement.regions.end(); ++it) {
//		BoundaryRegionStore *region = _mesh.bregion(it->first);
//		it->second.evaluator->evalSelected(
//				region->uniqueNodes->datatarray().size(),
//				region->uniqueNodes->datatarray().data(),
//				3, reinterpret_cast<double*>(_mesh.nodes->coordinates->datatarray().data()),
//				NULL, _step.currentTime, values.data() + offset);
//		offset += region->uniqueNodes->datatarray().size();
//	}
//
//	for (auto it = _stepSettings.displacement.intersections.begin(); it != _stepSettings.displacement.intersections.end(); ++it) {
//		BoundaryRegionsIntersectionStore *region = _mesh.ibregion(it->first);
//		it->second.evaluator->evalSelected(
//				region->uniqueNodes->datatarray().size(),
//				region->uniqueNodes->datatarray().data(),
//				3, reinterpret_cast<double*>(_mesh.nodes->coordinates->datatarray().data()),
//				NULL, _step.currentTime, values.data() + offset);
//		offset += region->uniqueNodes->datatarray().size();
//	}
}

void StructuralMechanics3DControler::initData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto c = _ncoordinate.data->begin(t);
		for (auto n = _mesh.elements->procNodes->datatarray().begin(t); n != _mesh.elements->procNodes->datatarray().end(t); ++n, ++c) {
			c->at(0) = _mesh.nodes->coordinates->datatarray()[*n].x;
			c->at(1) = _mesh.nodes->coordinates->datatarray()[*n].y;
			c->at(2) = _mesh.nodes->coordinates->datatarray()[*n].z;
		}
	}

	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = _step.currentTime;

	updateERegions(_globalSettings.initial_temperature, _nInitialTemperature.data->datatarray(), 1, cbegin, tbegin, time);

	for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = _mesh.boundaryRegions[r];
		if (region->dimension == 2) { // TODO: implement edge processing

			auto &distribution = region->procNodes->datatarray().distribution();

			_boundaries[r].coordinate.data = new serializededata<eslocal, double>(3, distribution);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto c = _boundaries[r].coordinate.data->begin(t);
				for (auto n = region->procNodes->datatarray().begin(t); n != region->procNodes->datatarray().end(t); ++n, ++c) {
					c->at(0) = _mesh.nodes->coordinates->datatarray()[*n].x;
					c->at(1) = _mesh.nodes->coordinates->datatarray()[*n].y;
					c->at(2) = _mesh.nodes->coordinates->datatarray()[*n].z;
				}
			}

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();

			auto pressure = _stepSettings.normal_pressure.find(region->name);
			if (pressure != _stepSettings.normal_pressure.end()) {
				updateBRegions(pressure->second, _boundaries[r].normalPressure, distribution, 2, cbegin, tbegin, time);
			}
		}
	}
}

void StructuralMechanics3DControler::nextTime()
{
	if (_step.isInitial()) {
		return;
	}

	parametersChanged();
}

void StructuralMechanics3DControler::parametersChanged()
{
	size_t threads = environment->OMP_NUM_THREADS;

	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = _step.currentTime;

	for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = _mesh.boundaryRegions[r];
		if (region->dimension == 2) {

			auto &distribution = region->procNodes->datatarray().distribution();

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();

			auto pressure = _stepSettings.normal_pressure.find(region->name);
			if (pressure != _stepSettings.normal_pressure.end()) {
				updateBRegions(pressure->second, _boundaries[r].normalPressure, distribution, 2, cbegin, tbegin, time);
			}
		}
	}
}

void StructuralMechanics3DControler::processElements(Matrices matrices, InstanceFiller &filler)
{
	auto enodes = _mesh.elements->procNodes->cbegin() + filler.begin;
	StructuralMechanics3DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - _mesh.elements->procNodes->datatarray().begin();
	iterator.temperature = _ntemperature.data->datatarray().begin() + noffset;
	iterator.coordinates = _ncoordinate.data->datatarray().begin() + noffset * 3;

	for (eslocal e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = _mesh.elements->epointers->datatarray()[e];
		iterator.material = _mesh.materials[_mesh.elements->material->datatarray()[e]];

		_kernel->processElement(matrices, iterator, _step, filler.Ke, filler.Me, filler.Re, filler.fe);
		filler.insert(enodes->size());

		iterator.temperature += enodes->size();
		iterator.coordinates += enodes->size() * 3;
	}
}

void StructuralMechanics3DControler::processBoundary(Matrices matrices, size_t rindex, InstanceFiller &filler)
{
	if (_mesh.boundaryRegions[rindex]->dimension != 1) {
		return;
	}

	auto enodes = _mesh.boundaryRegions[rindex]->procNodes->cbegin() + filler.begin;
	StructuralMechanics3DKernel::BoundaryIterator iterator;

	size_t noffset = enodes->begin() - _mesh.boundaryRegions[rindex]->procNodes->datatarray().begin();
	iterator.coordinates = _boundaries[rindex].coordinate.data->datatarray().begin() + noffset * 3;

	iterator.normalPressure = _boundaries[rindex].normalPressure.data ? _boundaries[rindex].normalPressure.data->datatarray().begin() + noffset : NULL;

	for (eslocal e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = _mesh.boundaryRegions[rindex]->epointers->datatarray()[e];

		_kernel->processEdge(matrices, iterator, _step, filler.Ke, filler.fe);
		filler.insert(enodes->size());

		if (iterator.normalPressure) {
			iterator.normalPressure += enodes->size();
		}
	}
}

void StructuralMechanics3DControler::processSolution()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		auto enodes = _mesh.elements->procNodes->cbegin(t);
		StructuralMechanics3DKernel::SolutionIterator iterator;

		size_t noffset = enodes->begin() - _mesh.elements->procNodes->datatarray().begin(t);
		iterator.temperature = _ntemperature.data->datatarray().begin(t);
		iterator.coordinates = _ncoordinate.data->datatarray().begin(t);

		for (size_t e = _mesh.elements->distribution[t]; e < _mesh.elements->distribution[t + 1]; ++e, ++enodes) {
			iterator.element = _mesh.elements->epointers->datatarray()[e];
			iterator.material = _mesh.materials[_mesh.elements->material->datatarray()[e]];

			_kernel->processSolution(iterator, _step);

			iterator.temperature += enodes->size();
			iterator.coordinates += enodes->size() * 3;
		}
	}
}







