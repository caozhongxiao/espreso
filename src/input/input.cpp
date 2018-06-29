
#include "input.h"
#include "workbench/workbench.h"
#include "../old/input/loader.h"

#include "../basis/containers/serializededata.h"
#include "../basis/logging/logging.h"
#include "../basis/utilities/communication.h"
#include "../basis/utilities/utils.h"

#include "../config/ecf/root.h"
#include "../mesh/mesh.h"
#include "../mesh/elements/element.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/elementsregionstore.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void Input::load(const ECFRoot &configuration, Mesh &mesh)
{
	switch (configuration.input) {
	case INPUT_FORMAT::WORKBENCH:
		WorkbenchLoader::load(configuration, mesh);
		mesh.update();
		break;
	default:
		input::OldLoader::load(configuration, *mesh.mesh, environment->MPIrank, environment->MPIsize);
		mesh.load();
		break;
	}
}

void Input::balance()
{
	if (environment->MPIsize == 1) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: use the sequential input for building mesh on 1 MPI process.";
	}

	int sorted, allSorted;

	sorted = std::is_sorted(_meshData.nIDs.begin(), _meshData.nIDs.end());
	MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

	if (allSorted) {
		balanceNodes();
	} else {
		balancePermutedNodes();
	}

	sorted = std::is_sorted(_meshData.eIDs.begin(), _meshData.eIDs.end());
	MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

	if (allSorted) {
		balanceElements();
	} else {
		balancePermutedElements();
	}
}

void Input::balanceNodes()
{
	std::vector<size_t> cCurrent = Communication::getDistribution(_meshData.nIDs.size(), MPITools::operations().sizeToOffsetsSize_t);
	_nDistribution = tarray<eslocal>::distribute(environment->MPIsize, cCurrent.back());

	if (!Communication::balance(_meshData.nIDs, cCurrent, _nDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance node IDs.";
	}
	if (!Communication::balance(_meshData.coordinates, cCurrent, _nDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance coordinates.";
	}
}

void Input::balancePermutedNodes()
{
	// TODO: optimize all to all

	eslocal myMaxID = 0, maxID;
	std::vector<eslocal> permutation(_meshData.nIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _meshData.nIDs[i] < _meshData.nIDs[j]; });

	if (_meshData.nIDs.size()) {
		myMaxID = _meshData.nIDs[permutation.back()];
	}
	MPI_Allreduce(&myMaxID, &maxID, sizeof(eslocal), MPI_BYTE, MPITools::operations().max, environment->MPICommunicator);

	_nDistribution = tarray<eslocal>::distribute(environment->MPIsize, maxID + 1);
	std::vector<std::vector<eslocal> > sIDs, rIDs;
	std::vector<std::vector<Point> > sCoordinates, rCoordinates;
	std::vector<int> targets;
	for (int r = 0; r < environment->MPIsize; r++) {
		auto begin = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r], [&] (eslocal i, const size_t &ID) { return _meshData.nIDs[i] < ID; });
		auto end = std::lower_bound(permutation.begin(), permutation.end(), _nDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _meshData.nIDs[i] < ID; });
		if (begin != end) {
			sIDs.push_back({});
			sCoordinates.push_back({});
			targets.push_back(r);
		}
		for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
			sIDs.back().push_back(_meshData.nIDs[permutation[n]]);
			sCoordinates.back().push_back(_meshData.coordinates[permutation[n]]);
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
	_meshData.nIDs.clear();
	_meshData.coordinates.clear();
	_meshData.nIDs.reserve(permutation.size());
	_meshData.coordinates.reserve(permutation.size());
	for (size_t n = 0; n < permutation.size(); n++) {
		_meshData.nIDs.push_back(rIDs[0][permutation[n]]);
		_meshData.coordinates.push_back(rCoordinates[0][permutation[n]]);
	}
}

void Input::balanceElements()
{
	std::vector<size_t> eCurrent = Communication::getDistribution(_meshData.esize.size(), MPITools::operations().sizeToOffsetsSize_t);
	std::vector<size_t> eTarget = tarray<eslocal>::distribute(environment->MPIsize, eCurrent.back());

	std::vector<size_t> nCurrent = Communication::getDistribution(_meshData.enodes.size(), MPITools::operations().sizeToOffsetsSize_t);
	std::vector<size_t> nTarget;

	if (environment->MPIrank == 0) {
		nTarget.push_back(0);
	}

	size_t nodeOffset = nCurrent[environment->MPIrank];
	size_t eTargetIndex = std::lower_bound(eTarget.begin(), eTarget.end(), eCurrent[environment->MPIrank] + 1) - eTarget.begin();
	for (size_t n = 0; n < _meshData.esize.size(); ++n) {
		nodeOffset += _meshData.esize[n];
		if (eCurrent[environment->MPIrank] + n + 1 == eTarget[eTargetIndex]) {
			nTarget.push_back(nodeOffset);
			++eTargetIndex;
		}
	}
	Communication::allGatherUnknownSize(nTarget);
	nTarget.resize(environment->MPIsize + 1, nTarget.back());

	if (!Communication::balance(_meshData.enodes, nCurrent, nTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element nodes.";
	}
	if (!Communication::balance(_meshData.esize, eCurrent, eTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element sizes.";
	}
	if (!Communication::balance(_meshData.eIDs, eCurrent, eTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element IDs.";
	}
	if (!Communication::balance(_meshData.body, eCurrent, eTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element bodies.";
	}
	if (!Communication::balance(_meshData.etype, eCurrent, eTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element types.";
	}
	if (!Communication::balance(_meshData.material, eCurrent, eTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element material.";
	}

	_eDistribution = eTarget;
}

void Input::balancePermutedElements()
{
	// TODO: optimize all to all

	eslocal myMaxID = 0, maxID;
	std::vector<eslocal> permutation(_meshData.esize.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _meshData.eIDs[i] < _meshData.eIDs[j]; });

	if (_meshData.esize.size()) {
		myMaxID = _meshData.eIDs[permutation.back()];
	}
	MPI_Allreduce(&myMaxID, &maxID, sizeof(eslocal), MPI_BYTE, MPITools::operations().max, environment->MPICommunicator);

	_eDistribution = tarray<eslocal>::distribute(environment->MPIsize, maxID + 1);
	std::vector<std::vector<eslocal> > sSize, sNodes, rSize, rNodes, sIDs, rIDs;
	std::vector<std::vector<int> > sType, rType, sMat, rMat, sBody, rBody;
	std::vector<int> targets;
	std::vector<eslocal> edist = { 0 };
	edist.reserve(_meshData.esize.size() + 1);
	for (size_t e = 0; e < _meshData.esize.size(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	for (int r = 0; r < environment->MPIsize; r++) {
		auto begin = std::lower_bound(permutation.begin(), permutation.end(), _eDistribution[r], [&] (eslocal i, const size_t &ID) { return _meshData.eIDs[i] < ID; });
		auto end = std::lower_bound(permutation.begin(), permutation.end(), _eDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _meshData.eIDs[i] < ID; });
		if (begin != end) {
			sSize.push_back({});
			sNodes.push_back({});
			sIDs.push_back({});
			sType.push_back({});
			sMat.push_back({});
			sBody.push_back({});
			targets.push_back(r);
		}
		for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
			sSize.back().push_back(_meshData.esize[permutation[n]]);
			sIDs.back().push_back(_meshData.eIDs[permutation[n]]);
			sType.back().push_back(_meshData.etype[permutation[n]]);
			sMat.back().push_back(_meshData.material[permutation[n]]);
			sBody.back().push_back(_meshData.body[permutation[n]]);
			sNodes.back().insert(sNodes.back().end(), _meshData.enodes.begin() + edist[permutation[n]], _meshData.enodes.begin() + edist[permutation[n] + 1]);
		}
	}

	if (!Communication::sendVariousTargets(sSize, rSize, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted elements sizes.";
	}
	if (!Communication::sendVariousTargets(sIDs, rIDs, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
	}
	if (!Communication::sendVariousTargets(sType, rType, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
	}
	if (!Communication::sendVariousTargets(sMat, rMat, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
	}
	if (!Communication::sendVariousTargets(sBody, rBody, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
	}
	if (!Communication::sendVariousTargets(sNodes, rNodes, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element nodes.";
	}

	for (size_t r = 1; r < rSize.size(); r++) {
		rSize[0].insert(rSize[0].end(), rSize[r].begin(), rSize[r].end());
		rIDs[0].insert(rIDs[0].end(), rIDs[r].begin(), rIDs[r].end());
		rType[0].insert(rType[0].end(), rType[r].begin(), rType[r].end());
		rMat[0].insert(rMat[0].end(), rMat[r].begin(), rMat[r].end());
		rBody[0].insert(rBody[0].end(), rBody[r].begin(), rBody[r].end());
		rNodes[0].insert(rNodes[0].end(), rNodes[r].begin(), rNodes[r].end());
	}

	permutation.resize(rSize[0].size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return rIDs[0][i] < rIDs[0][j]; });

	edist = std::vector<eslocal>({ 0 });
	edist.reserve(rSize[0].size() + 1);
	for (size_t e = 0; e < rSize[0].size(); e++) {
		edist.push_back(edist.back() + rSize[0][e]);
	}

	_meshData.esize.clear();
	_meshData.enodes.clear();
	_meshData.eIDs.clear();
	_meshData.etype.clear();
	_meshData.material.clear();
	_meshData.body.clear();
	_meshData.esize.reserve(permutation.size());
	_meshData.eIDs.reserve(permutation.size());
	_meshData.etype.reserve(permutation.size());
	_meshData.material.reserve(permutation.size());
	_meshData.body.reserve(permutation.size());
	_meshData.enodes.reserve(rNodes[0].size());
	for (size_t n = 0; n < permutation.size(); n++) {
		_meshData.esize.push_back(rSize[0][permutation[n]]);
		_meshData.eIDs.push_back(rIDs[0][permutation[n]]);
		_meshData.etype.push_back(rType[0][permutation[n]]);
		_meshData.material.push_back(rMat[0][permutation[n]]);
		_meshData.body.push_back(rBody[0][permutation[n]]);
		_meshData.enodes.insert(_meshData.enodes.end(), rNodes[0].begin() + edist[permutation[n]], rNodes[0].begin() + edist[permutation[n] + 1]);
	}
}

void Input::fillElements()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > tedist(threads);
	std::vector<std::vector<eslocal> > tnodes(threads);
	std::vector<std::vector<eslocal> > eIDs(threads), rData(threads);
	std::vector<std::vector<int> > eMat(threads), eBody(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	std::vector<eslocal> edist = { 0 };
	edist.reserve(_meshData.esize.size() + 1);
	for (size_t e = 0; e < _meshData.esize.size(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	std::vector<size_t> edistribution = tarray<Point>::distribute(threads, _meshData.esize.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if(t == 0) {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t], edist.begin() + edistribution[t + 1] + 1);
		} else {
			tedist[t].insert(tedist[t].end(), edist.begin() + edistribution[t] + 1, edist.begin() + edistribution[t + 1] + 1);
		}
		tnodes[t].insert(tnodes[t].end(), _meshData.enodes.begin() + edist[edistribution[t]], _meshData.enodes.begin() + edist[edistribution[t + 1]]);

		// till now, IDs are irelevant
		eIDs[t].resize(edistribution[t + 1] - edistribution[t]);
		std::iota(eIDs[t].begin(), eIDs[t].end(), _eDistribution[environment->MPIrank] + edistribution[t]);

		epointers[t].reserve(edistribution[t + 1] - edistribution[t]);
		eBody[t].reserve(edistribution[t + 1] - edistribution[t]);
		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
			epointers[t].push_back(&_mesh._eclasses[t][_meshData.etype[e]]);
			eBody[t].push_back(_meshData.body[e]);
		}

		if (_configuration.input == INPUT_FORMAT::WORKBENCH && _configuration.workbench.keep_material_sets) {
			eMat[t].reserve(edistribution[t + 1] - edistribution[t]);
			for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e) {
				eMat[t].push_back(_meshData.material[e]);
			}
		} else {
			eMat[t].resize(edistribution[t + 1] - edistribution[t]);
		}
	}

	_mesh.elements->dimension = _mesh.dimension;
	_mesh.elements->size = _meshData.esize.size();
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



