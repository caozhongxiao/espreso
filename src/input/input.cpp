
#include "input.h"
#include "workbench/workbench.h"
#include "../old/input/loader.h"

#include "../basis/containers/tarray.h"
#include "../basis/logging/logging.h"
#include "../basis/utilities/communication.h"
#include "../basis/utilities/utils.h"

#include "../config/ecf/root.h"
#include "../mesh/mesh.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void Input::load(const ECFRoot &configuration, Mesh &mesh, int MPIrank, int MPIsize)
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

void Input::balance()
{
	int sorted, allSorted;

	sorted = std::is_sorted(_meshData.nIDs.begin(), _meshData.nIDs.end());
	MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

	if (allSorted) {
		balanceNodes();
	} else {
		balancePermutedNodes();
	}

	sorted = std::is_sorted(_meshData.edata.begin(), _meshData.edata.end(), [] (const EData &e1, const EData &e2) { return e1.id < e2.id; });
	MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

	if (allSorted) {
		balanceElements();
	} else {
		balancePermutedElements();
	}
}

void Input::balanceNodes()
{
	if (environment->MPIsize == 1) {
		_nDistribution = { 0, _meshData.nIDs.size() };
		return;
	}

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

	if (environment->MPIsize == 1) {
		_nDistribution = { 0, _meshData.nIDs.size() };
		std::vector<eslocal> nIDs; nIDs.reserve(_meshData.nIDs.size());
		std::vector<Point> nPoints; nPoints.reserve(_meshData.nIDs.size());
		_meshData.nIDs.swap(nIDs);
		_meshData.coordinates.swap(nPoints);

		for (size_t i = 0; i < permutation.size(); i++) {
			_meshData.nIDs.push_back(nIDs[permutation[i]]);
			_meshData.coordinates.push_back(nPoints[permutation[i]]);
		}

		return;
	}

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
	if (environment->MPIsize == 1) {
		_eDistribution = { 0, _meshData.nIDs.size() };
		return;
	}

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

	if (!Communication::balance(_meshData.esize, eCurrent, eTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element sizes.";
	}
	if (!Communication::balance(_meshData.edata, eCurrent, eTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element data.";
	}
	if (!Communication::balance(_meshData.enodes, nCurrent, nTarget)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance element nodes.";
	}

	_eDistribution = eTarget;
}

void Input::balancePermutedElements()
{
	// TODO: optimize all to all

	eslocal myMaxID = 0, maxID;
	std::vector<eslocal> permutation(_meshData.edata.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _meshData.edata[i].id < _meshData.edata[j].id; });

	if (environment->MPIsize == 1) {
		_eDistribution = { 0, _meshData.esize.size() };
		std::vector<eslocal> eSize; eSize.reserve(_meshData.esize.size());
		std::vector<eslocal> eNodes; eNodes.reserve(_meshData.enodes.size());
		std::vector<EData> eData; eData.reserve(_meshData.esize.size());

		_meshData.esize.swap(eSize);
		_meshData.edata.swap(eData);
		_meshData.enodes.swap(eNodes);

		std::vector<eslocal> edist = { 0 };
		edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			edist.push_back(edist.back() + _meshData.esize[e]);
		}

		for (size_t i = 0; i < permutation.size(); i++) {
			_meshData.edata.push_back(eData[permutation[i]]);
			_meshData.esize.push_back(eSize[permutation[i]]);
			_meshData.enodes.insert(eNodes.end(), eNodes.begin() + edist[permutation[i]], eNodes.begin() + edist[permutation[i] + 1]);
		}

		return;
	}

	if (_meshData.esize.size()) {
		myMaxID = _meshData.edata[permutation.back()].id;
	}
	MPI_Allreduce(&myMaxID, &maxID, sizeof(eslocal), MPI_BYTE, MPITools::operations().max, environment->MPICommunicator);

	_eDistribution = tarray<eslocal>::distribute(environment->MPIsize, maxID + 1);
	std::vector<std::vector<eslocal> > sSize, sNodes, rSize, rNodes;
	std::vector<std::vector<EData> > sEData, rEData;
	std::vector<int> targets;
	std::vector<eslocal> edist = { 0 };
	edist.reserve(_meshData.esize.size() + 1);
	for (size_t e = 0; e < _meshData.esize.size(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	for (int r = 0; r < environment->MPIsize; r++) {
		auto begin = std::lower_bound(permutation.begin(), permutation.end(), _eDistribution[r], [&] (eslocal i, const size_t &ID) { return _meshData.edata[i].id < ID; });
		auto end = std::lower_bound(permutation.begin(), permutation.end(), _eDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _meshData.edata[i].id < ID; });
		if (begin != end) {
			sSize.push_back({});
			sNodes.push_back({});
			sEData.push_back({});
			targets.push_back(r);
		}
		for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
			sSize.back().push_back(_meshData.esize[permutation[n]]);
			sEData.back().push_back(_meshData.edata[permutation[n]]);
			sNodes.back().insert(sNodes.back().end(), _meshData.enodes.begin() + edist[permutation[n]], _meshData.enodes.begin() + edist[permutation[n] + 1]);
		}
	}

	if (!Communication::sendVariousTargets(sSize, rSize, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted elements sizes.";
	}
	if (!Communication::sendVariousTargets(sEData, rEData, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element data.";
	}
	if (!Communication::sendVariousTargets(sNodes, rNodes, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute not sorted element nodes.";
	}

	for (size_t r = 1; r < rSize.size(); r++) {
		rSize[0].insert(rSize[0].end(), rSize[r].begin(), rSize[r].end());
		rEData[0].insert(rEData[0].end(), rEData[r].begin(), rEData[r].end());
		rNodes[0].insert(rNodes[0].end(), rNodes[r].begin(), rNodes[r].end());
	}

	permutation.resize(rSize[0].size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return rEData[0][i].id < rEData[0][j].id; });

	edist = std::vector<eslocal>({ 0 });
	edist.reserve(rSize[0].size() + 1);
	for (size_t e = 0; e < rSize[0].size(); e++) {
		edist.push_back(edist.back() + rSize[0][e]);
	}

	_meshData.esize.clear();
	_meshData.enodes.clear();
	_meshData.edata.clear();
	_meshData.esize.reserve(permutation.size());
	_meshData.edata.reserve(permutation.size());
	_meshData.enodes.reserve(rNodes[0].size());
	for (size_t n = 0; n < permutation.size(); n++) {
		_meshData.esize.push_back(rSize[0][permutation[n]]);
		_meshData.edata.push_back(rEData[0][permutation[n]]);
		_meshData.enodes.insert(_meshData.enodes.end(), rNodes[0].begin() + edist[permutation[n]], rNodes[0].begin() + edist[permutation[n] + 1]);
	}
}



