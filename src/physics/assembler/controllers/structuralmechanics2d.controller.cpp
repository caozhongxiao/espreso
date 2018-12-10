
#include "structuralmechanics2d.controller.h"
#include "../kernels/structuralmechanics2d.kernel.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../config/ecf/environment.h"
#include "../../../config/ecf/physics/structuralmechanics.h"
#include "../../../globals/time.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/nodestore.h"
#include "../../dataholder.h"

using namespace espreso;

StructuralMechanics2DControler::StructuralMechanics2DControler(
		Mesh &mesh,
		const StructuralMechanicsGlobalSettings &gSettings,
		const StructuralMechanicsStepSettings &sSettings,
		const StructuralMechanicsOutputSettings &oSettings)
: StructuralMechanicsControler(mesh, gSettings, sSettings, oSettings)
{
	_kernel = new StructuralMechanics2DKernel(_globalSettings, _outputSettings);

	Point defaultMotion(0, 0, 0);
	double defaultHeat = 0; //273.15;
	double defaultThickness = 1;

	_ncoordinate.data = new serializededata<eslocal, double>(2, _nDistribution);
	_ntemperature.data = new serializededata<eslocal, double>(1, _nDistribution);
	_nInitialTemperature.data = new serializededata<eslocal, double>(1, _nDistribution);

	_nthickness.isConts = setDefault(gSettings.thickness, defaultThickness);
	_nthickness.data = new serializededata<eslocal, double>(1, _nDistribution, defaultThickness);

	_avgThickness = _mesh.nodes->appendData(1, { }); // printed on elements

	_boundaries.resize(_mesh.boundaryRegions.size());
}

StructuralMechanics2DControler::~StructuralMechanics2DControler()
{
	delete _kernel;
}

void StructuralMechanics2DControler::analyticRegularization(size_t domain, bool ortogonalCluster)
{
//	Point center; size_t np; double norm;
//	if (ortogonalCluster) {
//		center = _cCenter[_mesh->elements->clusters[domain]];
//		np = _cNp[_mesh->elements->clusters[domain]];
//		norm = _cNorm[_mesh->elements->clusters[domain]].x;
//	} else {
//		center = _dCenter[domain];
//		np = _dNp[domain];
//		norm = _dNorm[domain].x;
//	}
//
//	_instance->N1[domain].rows = _instance->domainDOFCount[domain];
//	_instance->N1[domain].cols = 3;
//	_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
//	_instance->N1[domain].type = 'G';
//
//	_instance->N1[domain].dense_values.reserve(_instance->N1[domain].nnz);
//
//	for (size_t c = 0; c < 2; c++) {
//		std::vector<double> kernel = { 0, 0 };
//		kernel[c] = 1 / std::sqrt(np);
//		for (size_t i = 0; i < _instance->domainDOFCount[domain] / 2; i++) {
//			_instance->N1[domain].dense_values.insert(_instance->N1[domain].dense_values.end(), kernel.begin(), kernel.end());
//		}
//	}
//
//	for (size_t i = 0; i < _mesh->nodes->dintervals[domain].size(); i++) {
//		for (eslocal n = _mesh->nodes->dintervals[domain][i].begin; n < _mesh->nodes->dintervals[domain][i].end; ++n) {
//			Point p = _mesh->nodes->coordinates->datatarray()[n] - center;
//			_instance->N1[domain].dense_values.push_back(-p.y / norm);
//			_instance->N1[domain].dense_values.push_back( p.x / norm);
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
//	Nt.rows = 3;
//	Nt.cols = _instance->K[domain].cols;
//	Nt.nnz  = 4 * fixPoints.size();
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
//	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
//
//	auto n2DOF = [&] (eslocal node) {
//		auto dit = _mesh->nodes->dintervals[domain].begin();
//		while (dit->end < node) { ++dit; }
//		return dit->DOFOffset + node - dit->begin;
//	};
//
//	for (size_t c = 0; c < 2; c++) {
//		for (size_t i = 0; i < fixPoints.size(); i++) {
//			COLS.push_back(2 * n2DOF(fixPoints[i]) + c + 1);
//		}
//	}
//	VALS.insert(VALS.end(), 2 * fixPoints.size(), 1);
//
//	for (size_t i = 0; i < fixPoints.size(); i++) {
//		const Point &p = _mesh->nodes->coordinates->datatarray()[fixPoints[i]];
//		COLS.push_back(2 * n2DOF(fixPoints[i]) + 0 + 1);
//		COLS.push_back(2 * n2DOF(fixPoints[i]) + 1 + 1);
//		VALS.push_back(-p.y);
//		VALS.push_back( p.x);
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

void StructuralMechanics2DControler::dirichletIndices(std::vector<std::vector<eslocal> > &indices)
{
	indices.resize(2);

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

void StructuralMechanics2DControler::dirichletValues(std::vector<double> &values)
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

void StructuralMechanics2DControler::initData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto c = _ncoordinate.data->begin(t);
		for (auto n = _mesh.elements->procNodes->datatarray().begin(t); n != _mesh.elements->procNodes->datatarray().end(t); ++n, ++c) {
			c->at(0) = _mesh.nodes->coordinates->datatarray()[*n].x;
			c->at(1) = _mesh.nodes->coordinates->datatarray()[*n].y;
		}
	}

	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = time::current;

	updateERegions(_globalSettings.initial_temperature, _nInitialTemperature.data->datatarray(), 1, cbegin, tbegin, time);
	updateERegions(_globalSettings.thickness, _nthickness.data->datatarray(), 1, cbegin, tbegin, time);

	averageNodeInitilization(_nthickness.data->datatarray(), _avgThickness->data);

	for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = _mesh.boundaryRegions[r];
		if (region->dimension == 1) {

			auto &distribution = region->procNodes->datatarray().distribution();

			_boundaries[r].coordinate.data = new serializededata<eslocal, double>(2, distribution);
			_boundaries[r].thickness.data = new serializededata<eslocal, double>(1, distribution);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto c = _boundaries[r].coordinate.data->begin(t);
				auto thick = _boundaries[r].thickness.data->begin(t);
				for (auto n = region->procNodes->datatarray().begin(t); n != region->procNodes->datatarray().end(t); ++n, ++c, ++thick) {
					c->at(0) = _mesh.nodes->coordinates->datatarray()[*n].x;
					c->at(1) = _mesh.nodes->coordinates->datatarray()[*n].y;
					thick->at(0) = _avgThickness->data[*n];
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

void StructuralMechanics2DControler::nextTime()
{
	if (time::isInitial()) {
		return;
	}

	parametersChanged();
}

void StructuralMechanics2DControler::parametersChanged()
{
	size_t threads = environment->OMP_NUM_THREADS;

	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = time::current;

	updateERegions(_globalSettings.thickness, _nthickness.data->datatarray(), 2, cbegin, tbegin, time);

	averageNodeInitilization(_nthickness.data->datatarray(), _avgThickness->data);

	for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = _mesh.boundaryRegions[r];
		if (region->dimension == 1) {

			auto &distribution = region->procNodes->datatarray().distribution();

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto thick = _boundaries[r].thickness.data->begin(t);
				for (auto n = region->procNodes->datatarray().begin(t); n != region->procNodes->datatarray().end(t); ++n, ++thick) {
					thick->at(0) = _avgThickness->data[*n];
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

void StructuralMechanics2DControler::processElements(Matrices matrices, InstanceFiller &filler)
{
	auto enodes = _mesh.elements->procNodes->cbegin() + filler.begin;
	StructuralMechanics2DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - _mesh.elements->procNodes->datatarray().begin();
	iterator.temperature = _ntemperature.data->datatarray().begin() + noffset;
	iterator.coordinates = _ncoordinate.data->datatarray().begin() + noffset * 2;
	iterator.thickness   = _nthickness.data->datatarray().begin() + noffset;

	for (eslocal e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = _mesh.elements->epointers->datatarray()[e];
		iterator.material = _mesh.materials[_mesh.elements->material->datatarray()[e]];

		_kernel->processElement(matrices, iterator, filler.Ke, filler.Me, filler.Re, filler.fe);
		filler.insert(enodes->size());

		iterator.temperature += enodes->size();
		iterator.coordinates += enodes->size() * 2;
		iterator.thickness   += enodes->size();
	}
}

void StructuralMechanics2DControler::processBoundary(Matrices matrices, size_t rindex, InstanceFiller &filler)
{
	if (_mesh.boundaryRegions[rindex]->dimension != 1) {
		return;
	}

	auto enodes = _mesh.boundaryRegions[rindex]->procNodes->cbegin() + filler.begin;
	StructuralMechanics2DKernel::BoundaryIterator iterator;

	size_t noffset = enodes->begin() - _mesh.boundaryRegions[rindex]->procNodes->datatarray().begin();
	iterator.coordinates = _boundaries[rindex].coordinate.data->datatarray().begin() + noffset * 2;
	iterator.thickness   = _boundaries[rindex].thickness.data->datatarray().begin() + noffset;

	iterator.normalPressure = _boundaries[rindex].normalPressure.data ? _boundaries[rindex].normalPressure.data->datatarray().begin() + noffset : NULL;

	for (eslocal e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = _mesh.boundaryRegions[rindex]->epointers->datatarray()[e];

		_kernel->processEdge(matrices, iterator, filler.Ke, filler.fe);
		filler.insert(enodes->size());

		iterator.coordinates += enodes->size() * 2;
		iterator.thickness   += enodes->size();
		if (iterator.normalPressure) {
			iterator.normalPressure += enodes->size();
		}
	}
}

void StructuralMechanics2DControler::processSolution()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		auto enodes = _mesh.elements->procNodes->cbegin(t);
		StructuralMechanics2DKernel::SolutionIterator iterator;

		size_t noffset = enodes->begin() - _mesh.elements->procNodes->datatarray().begin(t);
		iterator.temperature = _ntemperature.data->datatarray().begin(t);
		iterator.coordinates = _ncoordinate.data->datatarray().begin(t);
		iterator.thickness   = _nthickness.data->datatarray().begin(t);

		for (size_t e = _mesh.elements->distribution[t]; e < _mesh.elements->distribution[t + 1]; ++e, ++enodes) {
			iterator.element = _mesh.elements->epointers->datatarray()[e];
			iterator.material = _mesh.materials[_mesh.elements->material->datatarray()[e]];

			_kernel->processSolution(iterator);

			iterator.temperature += enodes->size();
			iterator.coordinates += enodes->size() * 2;
			iterator.thickness   += enodes->size();
		}
	}
}

