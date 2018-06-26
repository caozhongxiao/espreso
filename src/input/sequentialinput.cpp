
#include "sequentialinput.h"

#include "../basis/containers/serializededata.h"
#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"
#include "../basis/logging/timeeval.h"

#include "../config/ecf/root.h"

#include "../mesh/mesh.h"
#include "../mesh/preprocessing/meshpreprocessing.h"
#include "../mesh/elements/element.h"
#include "../mesh/store/nodestore.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/elementsregionstore.h"
#include "../mesh/store/boundaryregionstore.h"
#include "../old/input/loader.h"

#include <numeric>
#include <algorithm>
#include <fstream>

using namespace espreso;

void SequentialInput::buildMesh(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
{
	SequentialInput(configuration, meshData, mesh);
}

SequentialInput::SequentialInput(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
: Input(configuration, meshData, mesh)
{
	ESINFO(OVERVIEW) << "Build mesh.";
	TimeEval timing("Load sequential mesh");
	timing.totalTime.startWithBarrier();

	TimeEvent tnodes("fill nodes"); tnodes.start();
	_nDistribution = { 0, _meshData.nIDs.size() };
	fillNodes();
	tnodes.end(); timing.addEvent(tnodes);
	ESINFO(PROGRESS2) << "Sequential loader:: nodes filled.";

	TimeEvent telements("fill elements"); telements.start();
	_eDistribution = { 0, _meshData.esize.size() };
	fillElements();
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Sequential loader:: elements filled.";

	TimeEvent tnregions("fill nodes regions"); tnregions.start();
	fillNodeRegions();
	tnregions.end(); timing.addEvent(tnregions);
	ESINFO(PROGRESS2) << "Sequential loader:: nodes regions filled.";

	TimeEvent tbregions("fill boundary regions"); tbregions.start();
	fillBoundaryRegions();
	tbregions.end(); timing.addEvent(tbregions);
	ESINFO(PROGRESS2) << "Sequential loader:: boundary regions filled.";

	TimeEvent teregions("fill element regions"); teregions.start();
	fillElementRegions();
	teregions.end(); timing.addEvent(teregions);
	ESINFO(PROGRESS2) << "Sequential loader:: element regions filled.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void SequentialInput::fillNodes()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> cdistribution = tarray<Point>::distribute(threads, _meshData.coordinates.size());

	_mesh.nodes->size = _meshData.coordinates.size();
	_mesh.nodes->distribution = cdistribution;
	_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, tarray<eslocal>(cdistribution, _meshData.nIDs));
	_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, tarray<Point>(cdistribution, _meshData.coordinates));
	_mesh.nodes->ranks = new serializededata<eslocal, int>(1, tarray<int>(threads, _meshData.nIDs.size()));

	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, tarray<eslocal>(cdistribution, _meshData.nIDs));

	if (!std::is_sorted(_meshData.nIDs.begin(), _meshData.nIDs.end())) {
		std::vector<eslocal> permutation(_meshData.nIDs.size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _meshData.nIDs[i] < _meshData.nIDs[j]; });
		_mesh.nodes->permute(permutation);
	}

	_mesh.neighboursWithMe.push_back(environment->MPIrank);
}

void SequentialInput::fillNodeRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion) {
		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(nregion->first, _mesh._eclasses));
		_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, { threads, nregion->second });

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (auto n = _mesh.boundaryRegions.back()->nodes->begin(t)->begin(); n != _mesh.boundaryRegions.back()->nodes->end(t)->begin(); ++n) {
				*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
			}
		}
	}
}

void SequentialInput::fillBoundaryRegions()
{
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	for (size_t i = 0; i < _meshData.bregions.size(); i++) {
//		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_meshData.bregions[i].name, _mesh._eclasses));
//		_mesh.boundaryRegions.back()->distribution = tarray<Point>::distribute(threads, _meshData.bregions[i].esize.size());
//
//		std::vector<eslocal> edist = { 0 };
//		edist.reserve(_meshData.bregions[i].esize.size() + 1);
//		for (size_t e = 0; e < _meshData.bregions[i].esize.size(); e++) {
//			edist.push_back(edist.back() + _meshData.bregions[i].esize[e]);
//		}
//
//		std::vector<Element*> epointers;
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			for (size_t n = _mesh.boundaryRegions.back()->distribution[t]; n < _mesh.boundaryRegions.back()->distribution[t + 1]; ++n) {
//				epointers[n] = &_mesh._eclasses[t][_meshData.bregions[i].edata[n].etype];
//			}
//		}
//
//		switch (epointers.front()->type) {
//		case Element::TYPE::PLANE:
//			_mesh.boundaryRegions.back()->dimension = 2;
//			break;
//		case Element::TYPE::LINE:
//			_mesh.boundaryRegions.back()->dimension = 1;
//			break;
//		default:
//			ESINFO(ERROR) << "ESPRESO Workbench parser: invalid boundary region type. Have to be 3D plane or 2D line.";
//		}
//
//		std::vector<size_t> distribution = _mesh.boundaryRegions.back()->distribution;
//		for (size_t t = 0; t < threads; t++) {
//			++distribution[t];
//		}
//		_mesh.boundaryRegions.back()->elements = new serializededata<eslocal, eslocal>({ distribution, edist }, {threads, _meshData.bregions[i].enodes });
//		_mesh.boundaryRegions.back()->epointers = new serializededata<eslocal, Element*>(1, { threads, epointers });
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			for (auto n = _mesh.boundaryRegions.back()->elements->begin(t)->begin(); n != _mesh.boundaryRegions.back()->elements->end(t)->begin(); ++n) {
//				*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
//			}
//		}
//	}
}

void SequentialInput::fillElementRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion) {
		_mesh.elementsRegions.push_back(new ElementsRegionStore(eregion->first));
		_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, { threads, eregion->second });
	}
}



