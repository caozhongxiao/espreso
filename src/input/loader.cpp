
#include "loader.h"
#include "workbench/workbench.h"

#include "../basis/containers/point.h"
#include "../basis/containers/serializededata.h"
#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"

#include "../config/ecf/ecf.h"

#include "../mesh/mesh.h"
#include "../mesh/elements/element.h"
#include "../mesh/store/nodestore.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/elementsregionstore.h"
#include "../mesh/store/boundaryregionstore.h"
#include "../old/input/loader.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

void Loader::load(const ECFConfiguration &configuration, Mesh &mesh, int MPIrank, int MPIsize)
{
	switch (configuration.input) {
	case INPUT_FORMAT::WORKBENCH:
		WorkbenchLoader::load(configuration, mesh);
		mesh.update();
		break;
	default:
		input::OldLoader::load(configuration, *mesh.mesh, MPIrank, MPIsize);
		mesh.load();
		break;
	}
}

void Loader::loadDistributedMesh(DistributedMesh &dMesh, Mesh &mesh, bool shrinkIndices)
{
	Loader(dMesh, mesh, shrinkIndices);
}

Loader::Loader(DistributedMesh &dMesh, Mesh &mesh, bool shrinkIndices): _dMesh(dMesh), _mesh(mesh)
{
	fillMesh();
	addNodeRegions();
	addBoundaryRegions();

	if (shrinkIndices) {
		this->shrinkIndices();
	}
}

void Loader::fillMesh()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > nIDs(threads), nRanks(threads);
	std::vector<std::vector<Point> > tcoordinates(threads);

	std::vector<std::vector<eslocal> > tedist(threads);
	std::vector<std::vector<eslocal> > tnodes(threads);
	std::vector<std::vector<eslocal> > eIDs(threads), eMat(threads), eBody(threads), rData(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	std::vector<size_t> cdistribution = tarray<Point>::distribute(threads, _dMesh.coordinates.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		tcoordinates[t].insert(tcoordinates[t].end(), _dMesh.coordinates.begin() + cdistribution[t], _dMesh.coordinates.begin() + cdistribution[t + 1]);
		nIDs[t].insert(nIDs[t].end(), _dMesh.nIDs.begin() + cdistribution[t], _dMesh.nIDs.begin() + cdistribution[t + 1]);
		nRanks[t].resize(cdistribution[t + 1] - cdistribution[t]);
	}

	std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _dMesh.edist.size() - 1);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if(t == 0) {
			tedist[t].insert(tedist[t].end(), _dMesh.edist.begin() + edistribution[t], _dMesh.edist.begin() + edistribution[t + 1] + 1);
		} else {
			tedist[t].insert(tedist[t].end(), _dMesh.edist.begin() + edistribution[t] + 1, _dMesh.edist.begin() + edistribution[t + 1] + 1);
		}
		tnodes[t].insert(tnodes[t].end(), _dMesh.enodes.begin() + _dMesh.edist[edistribution[t]], _dMesh.enodes.begin() + _dMesh.edist[edistribution[t + 1]]);
		eIDs[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(eIDs[t].begin(), eIDs[t].end(), edistribution[t]);
		epointers[t].reserve(edistribution[t + 1] - edistribution[t]);

		eBody[t].reserve(edistribution[t + 1] - edistribution[t]);
		eMat[t].reserve(edistribution[t + 1] - edistribution[t]);
		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
			epointers[t].push_back(&_mesh._eclasses[t][_dMesh.edata[e].etype]);
			eBody[t].push_back(_dMesh.edata[e].body);
			eMat[t].push_back(_dMesh.edata[e].material);
		}
	}

	_mesh.nodes->size = _dMesh.coordinates.size();
	_mesh.nodes->distribution = cdistribution;
	_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, nIDs);
	_mesh.nodes->ranks = new serializededata<eslocal, int>(1, nRanks);
	_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, tcoordinates);

	_mesh.elements->size = _dMesh.edist.size() - 1;
	_mesh.elements->distribution = edistribution;
	_mesh.elements->IDs = new serializededata<eslocal, eslocal>(1, eIDs);
	_mesh.elements->nodes = new serializededata<eslocal, eslocal>(tedist, tnodes);
	_mesh.elements->epointers = new serializededata<eslocal, Element*>(1, epointers);
	_mesh.elements->material = new serializededata<eslocal, int>(1, eMat);
	_mesh.elements->body = new serializededata<eslocal, int>(1, eBody);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		rData[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(rData[t].begin(), rData[t].end(), edistribution[t]);
	}
	_mesh.elementsRegions.push_back(new ElementsRegionStore("ALL_ELEMENTS"));
	_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, rData);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		rData[t].resize(cdistribution[t + 1] - cdistribution[t]);
		std::iota(rData[t].begin(), rData[t].end(), cdistribution[t]);
	}
	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rData);

	_mesh.neighboursWithMe.push_back(environment->MPIrank);
	std::sort(_mesh.neighboursWithMe.begin(), _mesh.neighboursWithMe.end());
}

void Loader::addNodeRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	for (size_t i = 0; i < _dMesh.nregions.size(); i++) {
		std::vector<std::vector<eslocal> > tnodes(threads);
		std::sort(_dMesh.nregions[i].nodes.begin(), _dMesh.nregions[i].nodes.end());
		std::vector<size_t> tdistribution = tarray<eslocal>::distribute(threads, _dMesh.nregions[i].nodes.size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			tnodes[t].insert(tnodes[t].end(), _dMesh.nregions[i].nodes.begin() + tdistribution[t], _dMesh.nregions[i].nodes.begin() + tdistribution[t + 1]);
		}

		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.nregions[i].name, _mesh._eclasses));
		_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, tnodes);
	}
}

void Loader::addBoundaryRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;


	for (size_t i = 0; i < _dMesh.bregions.size(); i++) {

		std::vector<std::vector<eslocal> > tedist(threads), tnodes(threads);
		std::vector<std::vector<Element*> > epointers(threads);
		std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _dMesh.bregions[i].edist.size() - 1);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			if(t == 0) {
				tedist[t].insert(tedist[t].end(), _dMesh.bregions[i].edist.begin() + edistribution[t], _dMesh.bregions[i].edist.begin() + edistribution[t + 1] + 1);
			} else {
				tedist[t].insert(tedist[t].end(), _dMesh.bregions[i].edist.begin() + edistribution[t] + 1, _dMesh.bregions[i].edist.begin() + edistribution[t + 1] + 1);
			}
			tnodes[t].insert(tnodes[t].end(), _dMesh.bregions[i].enodes.begin() + _dMesh.bregions[i].edist[edistribution[t]], _dMesh.bregions[i].enodes.begin() + _dMesh.bregions[i].edist[edistribution[t + 1]]);
			epointers[t].reserve(edistribution[t + 1] - edistribution[t]);

			for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
				epointers[t].push_back(&_mesh._eclasses[t][_dMesh.bregions[i].etypes[e]]);
			}
		}

		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.bregions[i].name, _mesh._eclasses));
		_mesh.boundaryRegions.back()->distribution = edistribution;
		_mesh.boundaryRegions.back()->elements = new serializededata<eslocal, eslocal>(tedist, tnodes);
		_mesh.boundaryRegions.back()->epointers = new serializededata<eslocal, Element*>(1, epointers);
		_mesh.boundaryRegions.back()->dimension = 2;
	}
}

void Loader::shrinkIndices()
{
	std::vector<eslocal> permutation(_mesh.nodes->IDs->datatarray().size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _mesh.nodes->IDs->datatarray()[i] < _mesh.nodes->IDs->datatarray()[j]; });
	_mesh.nodes->permute(permutation);

	std::vector<eslocal> shrink(_mesh.nodes->IDs->datatarray().back(), -1);

	for (auto n = _mesh.nodes->IDs->datatarray().begin(); n != _mesh.nodes->IDs->datatarray().end(); ++n) {
		shrink[*n] = n - _mesh.nodes->IDs->datatarray().begin();
	}

	for (auto n = _mesh.elements->nodes->datatarray().begin(); n != _mesh.elements->nodes->datatarray().end(); ++n) {
		*n = shrink[*n];
	}

	for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
		if (_mesh.boundaryRegions[r]->dimension > 1) {
			for (auto n = _mesh.boundaryRegions[r]->elements->datatarray().begin(); n != _mesh.boundaryRegions[r]->elements->datatarray().end(); ++n) {
				*n = shrink[*n];
			}
		} else {
			for (auto n = _mesh.boundaryRegions[r]->nodes->datatarray().begin(); n != _mesh.boundaryRegions[r]->nodes->datatarray().end(); ++n) {
				*n = shrink[*n];
			}
		}
	}
}
