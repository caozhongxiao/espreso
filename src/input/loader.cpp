
#include "loader.h"
#include "workbench/workbench.h"

#include "../basis/containers/point.h"
#include "../basis/containers/serializededata.h"
#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"

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
#include "../config/ecf/root.h"

using namespace espreso;

void Loader::load(const ECFRoot &configuration, Mesh &mesh, int MPIrank, int MPIsize)
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

void Loader::loadDistributedMesh(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh)
{
	Loader(configuration, dMesh, mesh);
}

Loader::Loader(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh)
: _configuration(configuration), _dMesh(dMesh), _mesh(mesh)
{
	ESINFO(OVERVIEW) << "Balance distributed mesh.";
	distributeMesh();
	ESINFO(PROGRESS2) << "Distributed loader:: data balanced.";
	fillElements();
	ESINFO(PROGRESS2) << "Distributed loader:: elements filled.";
	fillCoordinates();
	ESINFO(PROGRESS2) << "Distributed loader:: coordinates filled.";
	addNodeRegions();
	ESINFO(PROGRESS2) << "Distributed loader:: node regions filled.";
	addBoundaryRegions();
	ESINFO(PROGRESS2) << "Distributed loader:: boundary regions filled.";
	addElementRegions();
	ESINFO(PROGRESS2) << "Distributed loader:: elements regions filled.";
}

void Loader::distributeMesh()
{
	eslocal myMaxID = 0, maxID;
	int sorted = std::is_sorted(_dMesh.nIDs.begin(), _dMesh.nIDs.end()), allSorted;
	std::vector<eslocal> permutation;
	if (_dMesh.nIDs.size()) {
		if (sorted) {
			myMaxID = _dMesh.nIDs.back();
		} else {
			permutation.resize(_dMesh.nIDs.size());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.nIDs[i] < _dMesh.nIDs[j]; });
			myMaxID = _dMesh.nIDs[permutation.back()];
		}
	}

	MPI_Allreduce(&myMaxID, &maxID, sizeof(eslocal), MPI_BYTE, MPITools::operations().max, environment->MPICommunicator);
	MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

	if (allSorted) {
		if (environment->MPIsize == 1) {
			_nDistribution = { 0, _dMesh.nIDs.size() };
			_eDistribution = { 0, _dMesh.esize.size() };
			return;
		}

		std::vector<size_t> cCurrent = Communication::getDistribution(_dMesh.nIDs.size());
		_nDistribution = tarray<eslocal>::distribute(environment->MPIsize, cCurrent.back());

		if (!Communication::balance(_dMesh.nIDs, cCurrent, _nDistribution)) {
			ESINFO(ERROR) << "ESPRESO internal error: balance node IDs.";
		}
		if (!Communication::balance(_dMesh.coordinates, cCurrent, _nDistribution)) {
			ESINFO(ERROR) << "ESPRESO internal error: balance coordinates.";
		}
	} else {
		if (environment->MPIsize == 1) {
			_nDistribution = { 0, _dMesh.nIDs.size() };
			_eDistribution = { 0, _dMesh.esize.size() };
			std::vector<eslocal> nIDs;
			std::vector<Point> nPoints;
			nIDs.reserve(_dMesh.nIDs.size());
			nPoints.reserve(_dMesh.nIDs.size());
			for (size_t i = 0; i < permutation.size(); i++) {
				nIDs.push_back(_dMesh.nIDs[permutation[i]]);
				nPoints.push_back(_dMesh.coordinates[permutation[i]]);
			}
			_dMesh.nIDs.swap(nIDs);
			_dMesh.coordinates.swap(nPoints);
			return;
		}

		_nDistribution = tarray<eslocal>::distribute(environment->MPIsize, maxID + 1);
		std::vector<std::vector<eslocal> > sIDs, rIDs;
		std::vector<std::vector<Point> > sCoordinates, rCoordinates;
		std::vector<int> targets;
		for (int r = 0; r < environment->MPIsize; r++) {
			auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r], [&] (eslocal i, const size_t &ID) { return _dMesh.nIDs[i] < ID; });
			auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _dMesh.nIDs[i] < ID; });
			if (begin != end) {
				sIDs.push_back({});
				sCoordinates.push_back({});
				targets.push_back(r);
			}
			for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
				sIDs.back().push_back(_dMesh.nIDs[permutation[n]]);
				sCoordinates.back().push_back(_dMesh.coordinates[permutation[n]]);
			}
		}

		if (!Communication::sendVariousTargets(sIDs, rIDs, targets)) {
			ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted node IDs.";
		}
		if (!Communication::sendVariousTargets(sCoordinates, rCoordinates, targets)) {
			ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted node coordinates.";
		}

		for (size_t r = 1; r < rIDs.size(); r++) {
			rIDs[0].insert(rIDs[0].end(), rIDs[r].begin(), rIDs[r].end());
			rCoordinates[0].insert(rCoordinates[0].end(), rCoordinates[r].begin(), rCoordinates[r].end());
		}

		permutation.resize(rIDs[0].size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return rIDs[0][i] < rIDs[0][j]; });
		_dMesh.nIDs.clear();
		_dMesh.coordinates.clear();
		_dMesh.nIDs.reserve(permutation.size());
		_dMesh.coordinates.reserve(permutation.size());
		for (size_t n = 0; n < permutation.size(); n++) {
			_dMesh.nIDs.push_back(rIDs[0][permutation[n]]);
			_dMesh.coordinates.push_back(rCoordinates[0][permutation[n]]);
		}
	}

	std::vector<size_t> eCurrent = Communication::getDistribution(_dMesh.esize.size());
	std::vector<size_t> eTarget = tarray<eslocal>::distribute(environment->MPIsize, eCurrent.back());

	std::vector<size_t> nCurrent = Communication::getDistribution(_dMesh.enodes.size());
	std::vector<size_t> nTarget;

	if (environment->MPIrank == 0) {
		nTarget.push_back(0);
	}

	size_t nodeOffset = nCurrent[environment->MPIrank];
	size_t eTargetIndex = std::lower_bound(eTarget.begin(), eTarget.end(), eCurrent[environment->MPIrank] + 1) - eTarget.begin();
	for (size_t n = 0; n < _dMesh.esize.size(); ++n) {
		nodeOffset += _dMesh.esize[n];
		if (eCurrent[environment->MPIrank] + n + 1 == eTarget[eTargetIndex]) {
			nTarget.push_back(nodeOffset);
			++eTargetIndex;
		}
	}
	Communication::allGatherUnknownSize(nTarget);

	if (!Communication::balance(_dMesh.esize, eCurrent, eTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element sizes.";
	}
	if (!Communication::balance(_dMesh.edata, eCurrent, eTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element data.";
	}
	if (!Communication::balance(_dMesh.enodes, nCurrent, nTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element nodes.";
	}

	_eDistribution = eTarget;
}

void Loader::fillElements()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > tedist(threads);
	std::vector<std::vector<eslocal> > tnodes(threads);
	std::vector<std::vector<eslocal> > eIDs(threads), eMat(threads), eBody(threads), rData(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	std::vector<eslocal> edist = { 0 };
	edist.reserve(_dMesh.esize.size() + 1);
	for (size_t e = 0; e < _dMesh.esize.size(); e++) {
		edist.push_back(edist.back() + _dMesh.esize[e]);
	}

	std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _dMesh.esize.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if(t == 0) {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t], edist.begin() + edistribution[t + 1] + 1);
		} else {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t] + 1, edist.begin() + edistribution[t + 1] + 1);
		}
		tnodes[t].insert(tnodes[t].end(), _dMesh.enodes.begin() + edist[edistribution[t]], _dMesh.enodes.begin() + edist[edistribution[t + 1]]);
		eIDs[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(eIDs[t].begin(), eIDs[t].end(), _eDistribution[environment->MPIrank] + edistribution[t]);
		epointers[t].reserve(edistribution[t + 1] - edistribution[t]);

		eBody[t].reserve(edistribution[t + 1] - edistribution[t]);

		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
			epointers[t].push_back(&_mesh._eclasses[t][_dMesh.edata[e].etype]);
			eBody[t].push_back(_dMesh.edata[e].body);
		}


		if (_configuration.input == INPUT_FORMAT::WORKBENCH && _configuration.workbench.keep_material_sets) {
			eMat[t].reserve(edistribution[t + 1] - edistribution[t]);
			for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
				eMat[t].push_back(_dMesh.edata[e].material);
			}
		} else {
			eMat[t].resize(edistribution[t + 1] - edistribution[t]);
		}
	}

	_mesh.elements->size = _dMesh.esize.size();
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
}

void Loader::fillCoordinates()
{
	size_t threads = environment->OMP_NUM_THREADS;

	if (environment->MPIsize == 1) {
		std::vector<std::vector<Point> > tcoordinates(threads);
		std::vector<std::vector<eslocal> > nIDs(threads), rData(threads);

		std::vector<size_t> cdistribution = tarray<Point>::distribute(threads, _dMesh.coordinates.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			tcoordinates[t].insert(tcoordinates[t].end(), _dMesh.coordinates.begin() + cdistribution[t], _dMesh.coordinates.begin() + cdistribution[t + 1]);
			nIDs[t].insert(nIDs[t].end(), _dMesh.nIDs.begin() + cdistribution[t], _dMesh.nIDs.begin() + cdistribution[t + 1]);
		}

		_mesh.nodes->size = _dMesh.coordinates.size();
		_mesh.nodes->distribution = cdistribution;
		_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, nIDs);
		_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, tcoordinates);
		_mesh.nodes->ranks = new serializededata<eslocal, int>(1, tarray<int>(threads, _nDistribution.back()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			rData[t].resize(cdistribution[t + 1] - cdistribution[t]);
			std::iota(rData[t].begin(), rData[t].end(), cdistribution[t]);
		}
		_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
		_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rData);

		_mesh.neighboursWithMe.push_back(environment->MPIrank);
		return;
	}

	std::vector<eslocal> nodes(_mesh.elements->nodes->datatarray().begin(), _mesh.elements->nodes->datatarray().end());
	Esutils::sortAndRemoveDuplicity(nodes);

	std::vector<std::vector<eslocal> > sBuffer;
	std::vector<int> sRanks;

	for (int t = 0; t < environment->MPIsize; t++) {
		auto begin = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[t]);
		auto end = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[t + 1]);
		if (end - begin) {
			sBuffer.push_back(std::vector<eslocal>(begin, end));
			sRanks.push_back(t);
		}
	}

	if (!Communication::sendVariousTargets(sBuffer, _rankNodeMap, sRanks, _targetRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange neighbors.";
	}

	std::vector<size_t> ndistribution = tarray<Point>::distribute(threads, _dMesh.coordinates.size());
	std::vector<std::vector<std::vector<eslocal> > > backedData(threads, std::vector<std::vector<eslocal> >(_targetRanks.size()));
	std::vector<std::vector<Point> > backedCoordinates(_targetRanks.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> ranks, ranksOffset;
		std::vector<std::vector<eslocal>::const_iterator> rPointer(_targetRanks.size());
		for (size_t r = 0; r < _targetRanks.size(); r++) {
			rPointer[r] = std::lower_bound(_rankNodeMap[r].begin(), _rankNodeMap[r].end(), _nDistribution[environment->MPIrank] + ndistribution[t]);
		}
		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			ranks.clear();
			ranksOffset.clear();
			for (size_t r = 0; r < _targetRanks.size(); r++) {
				if (rPointer[r] != _rankNodeMap[r].end() && *rPointer[r] == _nDistribution[environment->MPIrank] + n) {
					ranksOffset.push_back(r);
					ranks.push_back(_targetRanks[r]);
					++rPointer[r];
				}
			}
			for (size_t r = 0; r < ranks.size(); r++) {
				backedData[t][ranksOffset[r]].push_back(ranksOffset.size());
				backedData[t][ranksOffset[r]].insert(backedData[t][ranksOffset[r]].end(), ranks.begin(), ranks.end());
			}
		}
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t r = 0; r < _targetRanks.size(); r++) {
			backedData[0][r].insert(backedData[0][r].end(), backedData[t][r].begin(), backedData[t][r].end());
		}
	}

	for (size_t r = 0; r < _targetRanks.size(); r++) {
		backedCoordinates[r].reserve(_rankNodeMap[r].size());
		for (size_t n = 0; n < _rankNodeMap[r].size(); ++n) {
			backedCoordinates[r].push_back(_dMesh.coordinates[_rankNodeMap[r][n] - _nDistribution[environment->MPIrank]]);
		}
	}

	std::vector<std::vector<eslocal> > nodeRanks(sRanks.size()), allnodes(threads);
	std::vector<std::vector<Point> > coordinates(sRanks.size());

	if (!Communication::sendVariousTargets(backedData[0], nodeRanks, _targetRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return node ranks.";
	}
	if (!Communication::sendVariousTargets(backedCoordinates, coordinates, _targetRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return coordinates.";
	}

	size_t csize = 0;
	for (size_t i = 0; i < coordinates.size(); i++) {
		csize += coordinates[i].size();
	}

	std::vector<size_t> distribution = tarray<Point>::distribute(threads, csize);
	std::vector<std::vector<eslocal> > tIDs(threads), rankDistribution(sRanks.size());
	std::vector<std::vector<int> > rankData(sRanks.size());
	rankDistribution.front().push_back(0);

	#pragma omp parallel for
	for (size_t r = 0; r < sRanks.size(); r++) {
		for (size_t n = 0; n < nodeRanks[r].size(); n += nodeRanks[r][n] + 1) {
			rankData[r].insert(rankData[r].end(), nodeRanks[r].begin() + n + 1, nodeRanks[r].begin() + n + 1 + nodeRanks[r][n]);
			rankDistribution[r].push_back(rankData[r].size());
		}
	}

	Esutils::threadDistributionToFullDistribution(rankDistribution);

	for (size_t i = threads; i < sRanks.size(); i++) {
		coordinates[threads - 1].insert(coordinates[threads - 1].end(), coordinates[i].begin(), coordinates[i].end());
		rankData[threads - 1].insert(rankData[threads - 1].end(), rankData[i].begin(), rankData[i].end());
		rankDistribution[threads - 1].insert(rankDistribution[threads - 1].end(), rankDistribution[i].begin(), rankDistribution[i].end());
	}
	for (size_t i = threads; i < sRanks.size(); i++) {
		sBuffer[threads - 1].insert(sBuffer[threads - 1].end(), sBuffer[i].begin(), sBuffer[i].end());
	}
	coordinates.resize(threads);
	sBuffer.resize(threads);
	rankData.resize(threads);
	rankDistribution.resize(threads);

	serializededata<eslocal, Point>::balance(1, coordinates, &distribution);
	serializededata<eslocal, eslocal>::balance(1, sBuffer, &distribution);
	serializededata<eslocal, eslocal>::balance(rankDistribution, rankData, &distribution);


	_mesh.nodes->size = distribution.back();
	_mesh.nodes->distribution = distribution;
	_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, sBuffer);
	_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, coordinates);
	_mesh.nodes->ranks = new serializededata<eslocal, int>(rankDistribution, rankData);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		allnodes[t].resize(distribution[t + 1] - distribution[t]);
		std::iota(allnodes[t].begin(), allnodes[t].end(), distribution[t]);
		Esutils::sortAndRemoveDuplicity(rankData[t]);
	}

	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, allnodes);

	for (size_t t = 0; t < threads; t++) {
		_mesh.neighboursWithMe.insert(_mesh.neighboursWithMe.end(), rankData[t].begin(), rankData[t].end());
	}
	Esutils::sortAndRemoveDuplicity(_mesh.neighboursWithMe);

	for (size_t n = 0; n < _mesh.neighboursWithMe.size(); n++) {
		if (_mesh.neighboursWithMe[n] != environment->MPIrank) {
			_mesh.neighbours.push_back(_mesh.neighboursWithMe[n]);
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = _mesh.elements->nodes->begin(t)->begin(); n != _mesh.elements->nodes->end(t)->begin(); ++n) {
			*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
		}
	}
}

void Loader::addNodeRegions()
{
	// assume sorted nodes !!
	size_t threads = environment->OMP_NUM_THREADS;

	if (environment->MPIsize == 1) {
		for (size_t i = 0; i < _dMesh.nregions.size(); i++) {
			_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.nregions[i].name, _mesh._eclasses));

			std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, _dMesh.nregions[i].nodes.size());
			std::vector<std::vector<eslocal> > tnodes(threads);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				tnodes[t].insert(tnodes[t].end(), _dMesh.nregions[i].nodes.begin() + distribution[t], _dMesh.nregions[i].nodes.begin() + distribution[t + 1]);
			}
			_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, tnodes);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (auto n = _mesh.boundaryRegions.back()->nodes->begin(t)->begin(); n != _mesh.boundaryRegions.back()->nodes->end(t)->begin(); ++n) {
					*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
				}
			}
		}
		return;
	}

	for (size_t i = 0; i < _dMesh.nregions.size(); i++) {
		std::sort(_dMesh.nregions[i].nodes.begin(), _dMesh.nregions[i].nodes.end());

		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
		std::vector<int> sRanks, tRanks;

		for (int t = 0; t < environment->MPIsize; t++) {
			auto begin = std::lower_bound(_dMesh.nregions[i].nodes.begin(), _dMesh.nregions[i].nodes.end(), _nDistribution[t]);
			auto end = std::lower_bound(_dMesh.nregions[i].nodes.begin(), _dMesh.nregions[i].nodes.end(), _nDistribution[t + 1]);
			if (end - begin) {
				sBuffer.push_back(std::vector<eslocal>(begin, end));
				sRanks.push_back(t);
			}
		}

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
		}

		sBuffer.clear();
		sBuffer.resize(_targetRanks.size());
		for (size_t r = 1; r < rBuffer.size(); r++) {
			rBuffer[0].insert(rBuffer[0].end(), rBuffer[r].begin(), rBuffer[r].end());
		}

		if (rBuffer.size()) {
			#pragma omp parallel for
			for (size_t t = 0; t < _targetRanks.size(); t++) {
				sBuffer[t].resize(rBuffer[0].size());
				sBuffer[t].resize(std::set_intersection(_rankNodeMap[t].begin(), _rankNodeMap[t].end(), rBuffer[0].begin(), rBuffer[0].end(), sBuffer[t].begin()) - sBuffer[t].begin());
			}
		}

		for (size_t t = 0; t < _targetRanks.size(); t++) {
			if (sBuffer[t].size()) {
				tRanks.push_back(t);
			}
		}
		for (size_t t = 0; t < tRanks.size(); t++) {
			sBuffer[t].swap(sBuffer[tRanks[t]]);
			tRanks[t] = _targetRanks[tRanks[t]];
		}
		sBuffer.resize(tRanks.size());

		rBuffer.clear();
		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region to targets.";
		}

		for (size_t t = threads; t < rBuffer.size(); t++) {
			rBuffer[threads - 1].insert(rBuffer[threads - 1].end(), rBuffer[t].begin(), rBuffer[t].end());
		}
		rBuffer.resize(threads);
		serializededata<eslocal, eslocal>::balance(1, rBuffer);

		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.nregions[i].name, _mesh._eclasses));
		_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rBuffer);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (auto n = _mesh.boundaryRegions.back()->nodes->begin(t)->begin(); n != _mesh.boundaryRegions.back()->nodes->end(t)->begin(); ++n) {
				*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
			}
		}
	}
}

void Loader::addBoundaryRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	if (environment->MPIsize == 1) {
		for (size_t i = 0; i < _dMesh.bregions.size(); i++) {
			std::vector<eslocal> edist = { 0 };
			edist.reserve(_dMesh.bregions[i].esize.size() + 1);
			for (size_t e = 0; e < _dMesh.bregions[i].esize.size(); e++) {
				edist.push_back(edist.back() + _dMesh.bregions[i].esize[e]);
			}

			std::vector<std::vector<eslocal> > tedist(threads), tnodes(threads);
			std::vector<std::vector<Element*> > epointers(threads);
			std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _dMesh.bregions[i].esize.size());

			tedist.front().push_back(0);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t n = edistribution[t]; n < edistribution[t + 1]; ++n) {
					tnodes[t].insert(tnodes[t].end(), _dMesh.bregions[i].enodes.begin() + edist[n], _dMesh.bregions[i].enodes.begin() + edist[n + 1]);
					epointers[t].push_back(&_mesh._eclasses[t][_dMesh.bregions[i].etypes[n]]);
					tedist[t].push_back(tnodes[t].size());
				}
			}

			Esutils::threadDistributionToFullDistribution(tedist);

			_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.bregions[i].name, _mesh._eclasses));
			_mesh.boundaryRegions.back()->distribution = tarray<eslocal>::distribute(threads, epointers.front().size());
			_mesh.boundaryRegions.back()->dimension = 2;
			_mesh.boundaryRegions.back()->elements = new serializededata<eslocal, eslocal>(tedist, tnodes);
			_mesh.boundaryRegions.back()->epointers = new serializededata<eslocal, Element*>(1, epointers);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (auto n = _mesh.boundaryRegions.back()->elements->begin(t)->begin(); n != _mesh.boundaryRegions.back()->elements->end(t)->begin(); ++n) {
					*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
				}
			}
		}
		return;
	}

	if (_dMesh.bregions.size()) {
		_mesh.preprocessing->linkNodesAndElements();
	}
	std::vector<eslocal> edistribution = _mesh.elements->gatherElementsProcDistribution();

	for (size_t i = 0; i < _dMesh.bregions.size(); i++) {

		std::vector<eslocal> edist = { 0 };
		edist.reserve(_dMesh.bregions[i].esize.size() + 1);
		for (size_t e = 0; e < _dMesh.bregions[i].esize.size(); e++) {
			edist.push_back(edist.back() + _dMesh.bregions[i].esize[e]);
		}

		std::vector<eslocal> permutation(edist.size() - 1);
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (eslocal e1, eslocal e2) {
			return _dMesh.bregions[i].enodes[edist[e1]] < _dMesh.bregions[i].enodes[edist[e2]];
		});


		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
		std::vector<int> sRanks, tRanks;

		for (int t = 0; t < environment->MPIsize; t++) {
			auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[t], [&] (eslocal e, eslocal n) { return _dMesh.bregions[i].enodes[edist[e]] < n; });
			auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[t + 1], [&] (eslocal e, eslocal n) { return _dMesh.bregions[i].enodes[edist[e]] < n; });
			if (begin != end) {
				sRanks.push_back(t);
			}
		}
		sBuffer.resize(sRanks.size());

		#pragma omp parallel for
		for (size_t t = 0; t < sRanks.size(); t++) {
			auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[sRanks[t]], [&] (eslocal e, eslocal n) { return _dMesh.bregions[i].enodes[edist[e]] < n; });
			auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[sRanks[t] + 1], [&] (eslocal e, eslocal n) { return _dMesh.bregions[i].enodes[edist[e]] < n; });
			for (auto e = begin; e != end; ++e) {
				sBuffer[t].push_back(_dMesh.bregions[i].etypes[*e]);
				sBuffer[t].push_back(_dMesh.bregions[i].esize[*e]);
				for (eslocal n = 0; n < _dMesh.bregions[i].esize[*e]; ++n) {
					sBuffer[t].push_back(_dMesh.bregions[i].enodes[edist[*e] + n]);
				}
			}
		}

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
		}

		sBuffer.clear();
		sBuffer.resize(_targetRanks.size());

		for (size_t r = 0; r < rBuffer.size(); r++) {
			std::vector<eslocal> nodes;
			for (size_t n = 0; n < rBuffer[r].size(); n += 2 + rBuffer[r][n + 1]) {
				nodes.clear();
				nodes.insert(nodes.end(), rBuffer[r].begin() + n + 2, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
				std::sort(nodes.begin(), nodes.end());
				auto nbegin = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[environment->MPIrank]);
				auto nend = std::lower_bound(nodes.begin(), nodes.end(), _nDistribution[environment->MPIrank + 1]);

				#pragma omp parallel for
				for (size_t t = 0; t < _targetRanks.size(); t++) {
					auto it = _rankNodeMap[t].begin();
					bool found = true;
					for (auto current = nbegin; found && current != nend; ++current) {
						it = std::lower_bound(it, _rankNodeMap[t].end(), *current);
						found = it != _rankNodeMap[t].end() && *it == *current;
					}
					if (found) {
						sBuffer[t].insert(sBuffer[t].end(), rBuffer[r].begin() + n, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
					}
				}
			}
		}

		for (size_t t = 0; t < _targetRanks.size(); t++) {
			if (sBuffer[t].size()) {
				tRanks.push_back(t);
			}
		}
		for (size_t t = 0; t < tRanks.size(); t++) {
			sBuffer[t].swap(sBuffer[tRanks[t]]);
			tRanks[t] = _targetRanks[tRanks[t]];
		}
		sBuffer.resize(tRanks.size());

		rBuffer.clear();
		if (!Communication::sendVariousTargets(sBuffer, rBuffer, tRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region to targets.";
		}

		std::vector<std::vector<eslocal> > tedist(std::max((size_t)1, rBuffer.size()), { 0 }), tnodes(std::max((size_t)1, rBuffer.size()));
		std::vector<std::vector<Element*> > epointers(std::max((size_t)1, rBuffer.size()));

		#pragma omp parallel for
		for (size_t r = 0; r < rBuffer.size(); r++) {
			std::vector<eslocal> nodes;
			for (size_t n = 0; n < rBuffer[r].size(); n += 2 + rBuffer[r][n + 1]) {
				nodes.clear();
				nodes.insert(nodes.end(), rBuffer[r].begin() + n + 2, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
				std::sort(nodes.begin(), nodes.end());

				bool found = true;
				auto it = _mesh.nodes->IDs->datatarray().begin();
				for (auto current = nodes.begin(); found && current != nodes.end(); ++current) {
					it = std::lower_bound(it, _mesh.nodes->IDs->datatarray().end(), *current);
					found = it != _mesh.nodes->IDs->datatarray().end() && *it == *current;
				}
				if (found) {
					tnodes[r].insert(tnodes[r].end(), rBuffer[r].begin() + n + 2, rBuffer[r].begin() + n + 2 + rBuffer[r][n + 1]);
					tedist[r].push_back(tnodes[r].size());
					epointers[r].push_back(&_mesh._eclasses[0][rBuffer[r][n]]);
				}
			}
		}

		#pragma omp parallel for
		for (size_t r = 0; r < rBuffer.size(); r++) {
			for (auto n = tnodes[r].begin(); n != tnodes[r].end(); ++n) {
				*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
			}

			auto elinks = _mesh.nodes->elements->begin();
			std::vector<eslocal> elements;
			eslocal size, maxsize, first = edistribution[environment->MPIrank], last = edistribution[environment->MPIrank + 1];
			eslocal pdist = r == 0 ? 1 : 0, pnodes = 0, pclass = 0;
			for (size_t e = 0; e < epointers[r].size(); ++e) {
				elements.clear();
				for (eslocal n = tedist[r][e]; n < tedist[r][e] + epointers[r][e]->coarseNodes; ++n) {
					elements.insert(elements.end(), (elinks + tnodes[r][n])->begin(), (elinks + tnodes[r][n])->end());
				}
				std::sort(elements.begin(), elements.end());
				size = 1; maxsize = 0;
				for (size_t ei = 1; ei < elements.size(); ++ei) {
					if (elements[ei - 1] == elements[ei]) {
						++size;
					} else {
						if (maxsize < size && first <= elements[ei - 1] && elements[ei - 1] < last) {
							maxsize = size;
						}
						size = 1;
					}
				}
				if (maxsize < size && first <= elements.back() && elements.back() < last) {
					maxsize = size;
				}
				if (maxsize == epointers[r][e]->coarseNodes) {
					for (eslocal n = tedist[r][e]; n < tedist[r][e + 1]; ++n) {
						tnodes[r][pnodes++] = tnodes[r][n];
					}
					tedist[r][pdist++] = pnodes;
					epointers[r][pclass++] = epointers[r][e];
				}
			}
			tnodes[r].resize(pnodes);
			tedist[r].resize(pdist);
			epointers[r].resize(pclass);
		}

		Esutils::threadDistributionToFullDistribution(tedist);

		for (size_t t = threads; t < rBuffer.size(); t++) {
			tedist[threads - 1].insert(tedist[threads - 1].end(), tedist[t].begin(), tedist[t].end());
			tnodes[threads - 1].insert(tnodes[threads - 1].end(), tnodes[t].begin(), tnodes[t].end());
			epointers[threads - 1].insert(epointers[threads - 1].end(), epointers[t].begin(), epointers[t].end());
		}

		tedist.resize(threads);
		tnodes.resize(threads);
		epointers.resize(threads);

		serializededata<eslocal, eslocal>::balance(tedist, tnodes);
		serializededata<eslocal, Element*>::balance(1, epointers);

		#pragma omp parallel for
		for (size_t t = 1; t < threads; t++) {
			for (size_t e = 0; e < epointers[t].size(); e++) {
				epointers[t][e] = &_mesh._eclasses[t][epointers[t][e] - _mesh._eclasses[0]];
			}
		}

		_mesh.boundaryRegions.push_back(new BoundaryRegionStore(_dMesh.bregions[i].name, _mesh._eclasses));
		_mesh.boundaryRegions.back()->dimension = 2;
		_mesh.boundaryRegions.back()->elements = new serializededata<eslocal, eslocal>(tedist, tnodes);
		_mesh.boundaryRegions.back()->epointers = new serializededata<eslocal, Element*>(1, epointers);
		_mesh.boundaryRegions.back()->distribution = _mesh.boundaryRegions.back()->epointers->datatarray().distribution();
	}
}

void Loader::addElementRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	if (environment->MPIsize == 1) {
		for (size_t i = 0; i < _dMesh.eregions.size(); i++) {
			_mesh.elementsRegions.push_back(new ElementsRegionStore(_dMesh.eregions[i].name));

			std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, _dMesh.eregions[i].elements.size());
			std::vector<std::vector<eslocal> > telements(threads);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				telements[t].insert(telements[t].end(), _dMesh.eregions[i].elements.begin() + distribution[t], _dMesh.eregions[i].elements.begin() + distribution[t + 1]);
			}
			_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, telements);
		}
		return;
	}

	for (size_t i = 0; i < _dMesh.eregions.size(); i++) {
		std::vector<std::vector<eslocal> > sBuffer, rBuffer;
		std::vector<int> sRanks, tRanks;

		for (int t = 0; t < environment->MPIsize; t++) {
			auto begin = std::lower_bound(_dMesh.eregions[i].elements.begin(), _dMesh.eregions[i].elements.end(), _eDistribution[t]);
			auto end = std::lower_bound(_dMesh.eregions[i].elements.begin(), _dMesh.eregions[i].elements.end(), _eDistribution[t + 1]);
			if (end - begin) {
				sBuffer.push_back(std::vector<eslocal>(begin, end));
				sRanks.push_back(t);
			}
		}

		if (!Communication::sendVariousTargets(sBuffer, rBuffer, sRanks)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange node region.";
		}

		for (size_t t = threads; t < rBuffer.size(); t++) {
			rBuffer[threads - 1].insert(rBuffer[threads - 1].end(), rBuffer[t].begin(), rBuffer[t].end());
		}
		rBuffer.resize(threads);
		serializededata<eslocal, eslocal>::balance(1, rBuffer);

		_mesh.elementsRegions.push_back(new ElementsRegionStore(_dMesh.eregions[i].name));
		_mesh.elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, rBuffer);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (auto e = _mesh.elementsRegions.back()->elements->begin(t)->begin(); e != _mesh.elementsRegions.back()->elements->end(t)->begin(); ++e) {
				*e -= _eDistribution[environment->MPIrank];
			}
		}
	}
}
