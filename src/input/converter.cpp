
#include "converter.h"
#include "workbench/workbench.h"

#include "../basis/containers/point.h"
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

void Converter::load(const ECFRoot &configuration, Mesh &mesh, int MPIrank, int MPIsize)
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

void Converter::loadDistributedMesh(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh)
{
	Converter(configuration, dMesh, mesh);
}

Converter::Converter(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh)
: _configuration(configuration), _dMesh(dMesh), _mesh(mesh), _sfc(_mesh.dimension, SFCDEPTH, _dMesh.coordinates)
{

	ESINFO(OVERVIEW) << "Balance distributed mesh.";
	TimeEval timing("Load distributed mesh");
	timing.totalTime.startWithBarrier();

	TimeEvent tdistribution("distribute mesh across processes"); tdistribution.start();
	balance();
	tdistribution.end(); timing.addEvent(tdistribution);
	ESINFO(PROGRESS2) << "Balanced loader:: data balanced.";

	TimeEvent tnbuckets("assign nodes to buckets"); tnbuckets.start();
	assignNBuckets();
	tnbuckets.end(); timing.addEvent(tnbuckets);
	ESINFO(PROGRESS2) << "Balanced loader:: nodes buckets assigned.";

	TimeEvent tebuckets("assign elements to buckets"); tebuckets.start();
	assignEBuckets();
	tebuckets.end(); timing.addEvent(tebuckets);
	ESINFO(PROGRESS2) << "Balanced loader:: elements buckets assigned.";

	TimeEvent tclusterization("clusterization"); tclusterization.start();
	clusterize();
	tclusterization.end(); timing.addEvent(tclusterization);
	ESINFO(PROGRESS2) << "Balanced loader:: elements clusterized.";

	TimeEvent tlinkup("link together"); tlinkup.start();
	linkup();
	tlinkup.end(); timing.addEvent(tlinkup);
	ESINFO(PROGRESS2) << "Balanced loader:: neighbors linked up.";

	TimeEvent telements("fill elements"); telements.start();
	fillElements();
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Balanced loader:: elements filled.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void Converter::balance()
{
	int sorted, allSorted;

	sorted = std::is_sorted(_dMesh.nIDs.begin(), _dMesh.nIDs.end());
	MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

	if (allSorted) {
		balanceNodes();
	} else {
		balancePermutedNodes();
	}

	sorted = std::is_sorted(_dMesh.edata.begin(), _dMesh.edata.end(), [] (const EData &e1, const EData &e2) { return e1.id < e2.id; });
	MPI_Allreduce(&sorted, &allSorted, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);

	if (allSorted) {
		balanceElements();
	} else {
		balancePermutedElements();
	}
}

void Converter::balanceNodes()
{
	if (environment->MPIsize == 1) {
		_nDistribution = { 0, _dMesh.nIDs.size() };
		return;
	}

	std::vector<size_t> cCurrent = Communication::getDistribution(_dMesh.nIDs.size(), MPITools::operations().sizeToOffsetsSize_t);
	_nDistribution = tarray<eslocal>::distribute(environment->MPIsize, cCurrent.back());

	if (!Communication::balance(_dMesh.nIDs, cCurrent, _nDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance node IDs.";
	}
	if (!Communication::balance(_dMesh.coordinates, cCurrent, _nDistribution)) {
		ESINFO(ERROR) << "ESPRESO internal error: balance coordinates.";
	}
}

void Converter::balancePermutedNodes()
{
	// TODO: optimize all to all

	eslocal myMaxID = 0, maxID;
	std::vector<eslocal> permutation(_dMesh.nIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.nIDs[i] < _dMesh.nIDs[j]; });

	if (environment->MPIsize == 1) {
		_nDistribution = { 0, _dMesh.nIDs.size() };
		std::vector<eslocal> nIDs; nIDs.reserve(_dMesh.nIDs.size());
		std::vector<Point> nPoints; nPoints.reserve(_dMesh.nIDs.size());
		_dMesh.nIDs.swap(nIDs);
		_dMesh.coordinates.swap(nPoints);

		for (size_t i = 0; i < permutation.size(); i++) {
			_dMesh.nIDs.push_back(nIDs[permutation[i]]);
			_dMesh.coordinates.push_back(nPoints[permutation[i]]);
		}

		return;
	}

	if (_dMesh.nIDs.size()) {
		myMaxID = _dMesh.nIDs[permutation.back()];
	}
	MPI_Allreduce(&myMaxID, &maxID, sizeof(eslocal), MPI_BYTE, MPITools::operations().max, environment->MPICommunicator);

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

void Converter::balanceElements()
{
	if (environment->MPIsize == 1) {
		_eDistribution = { 0, _dMesh.nIDs.size() };
		return;
	}

	std::vector<size_t> eCurrent = Communication::getDistribution(_dMesh.esize.size(), MPITools::operations().sizeToOffsetsSize_t);
	std::vector<size_t> eTarget = tarray<eslocal>::distribute(environment->MPIsize, eCurrent.back());

	std::vector<size_t> nCurrent = Communication::getDistribution(_dMesh.enodes.size(), MPITools::operations().sizeToOffsetsSize_t);
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
	nTarget.resize(environment->MPIsize + 1, nTarget.back());

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

void Converter::balancePermutedElements()
{
	// TODO: optimize all to all

	eslocal myMaxID = 0, maxID;
	std::vector<eslocal> permutation(_dMesh.edata.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.edata[i].id < _dMesh.edata[j].id; });

	if (environment->MPIsize == 1) {
		_eDistribution = { 0, _dMesh.esize.size() };
		std::vector<eslocal> eSize; eSize.reserve(_dMesh.esize.size());
		std::vector<eslocal> eNodes; eNodes.reserve(_dMesh.enodes.size());
		std::vector<EData> eData; eData.reserve(_dMesh.esize.size());

		_dMesh.esize.swap(eSize);
		_dMesh.edata.swap(eData);
		_dMesh.enodes.swap(eNodes);

		std::vector<eslocal> edist = { 0 };
		edist.reserve(_dMesh.esize.size() + 1);
		for (size_t e = 0; e < _dMesh.esize.size(); e++) {
			edist.push_back(edist.back() + _dMesh.esize[e]);
		}

		for (size_t i = 0; i < permutation.size(); i++) {
			_dMesh.edata.push_back(eData[permutation[i]]);
			_dMesh.esize.push_back(eSize[permutation[i]]);
			_dMesh.enodes.insert(eNodes.end(), eNodes.begin() + edist[permutation[i]], eNodes.begin() + edist[permutation[i] + 1]);
		}

		return;
	}

	if (_dMesh.esize.size()) {
		myMaxID = _dMesh.edata[permutation.back()].id;
	}
	MPI_Allreduce(&myMaxID, &maxID, sizeof(eslocal), MPI_BYTE, MPITools::operations().max, environment->MPICommunicator);

	_eDistribution = tarray<eslocal>::distribute(environment->MPIsize, maxID + 1);
	std::vector<std::vector<eslocal> > sSize, sNodes, rSize, rNodes;
	std::vector<std::vector<EData> > sEData, rEData;
	std::vector<int> targets;
	std::vector<eslocal> edist = { 0 };
	edist.reserve(_dMesh.esize.size() + 1);
	for (size_t e = 0; e < _dMesh.esize.size(); e++) {
		edist.push_back(edist.back() + _dMesh.esize[e]);
	}

	for (int r = 0; r < environment->MPIsize; r++) {
		auto begin = std::lower_bound(permutation.begin(), permutation.end(), _eDistribution[r], [&] (eslocal i, const size_t &ID) { return _dMesh.edata[i].id < ID; });
		auto end = std::lower_bound(permutation.begin(), permutation.end(), _eDistribution[r + 1], [&] (eslocal i, const size_t &ID) { return _dMesh.edata[i].id < ID; });
		if (begin != end) {
			sSize.push_back({});
			sNodes.push_back({});
			sEData.push_back({});
			targets.push_back(r);
		}
		for (size_t n = begin - permutation.begin(); n < end - permutation.begin(); ++n) {
			sSize.back().push_back(_dMesh.esize[permutation[n]]);
			sEData.back().push_back(_dMesh.edata[permutation[n]]);
			sNodes.back().insert(sNodes.back().end(), _dMesh.enodes.begin() + edist[permutation[n]], _dMesh.enodes.begin() + edist[permutation[n] + 1]);
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

	_dMesh.esize.clear();
	_dMesh.enodes.clear();
	_dMesh.edata.clear();
	_dMesh.esize.reserve(permutation.size());
	_dMesh.edata.reserve(permutation.size());
	_dMesh.enodes.reserve(rNodes[0].size());
	for (size_t n = 0; n < permutation.size(); n++) {
		_dMesh.esize.push_back(rSize[0][permutation[n]]);
		_dMesh.edata.push_back(rEData[0][permutation[n]]);
		_dMesh.enodes.insert(_dMesh.enodes.end(), rNodes[0].begin() + edist[permutation[n]], rNodes[0].begin() + edist[permutation[n] + 1]);
	}
}

void Converter::assignNBuckets()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> cdistribution = tarray<Point>::distribute(threads, _dMesh.coordinates.size());

	_nBuckets.resize(_dMesh.coordinates.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t n = cdistribution[t]; n < cdistribution[t + 1]; ++n) {
			_nBuckets[n] = _sfc.getBucket(_dMesh.coordinates[n]);
		}
	}
}

void Converter::assignEBuckets()
{ // we needs to ask a neighbor process to get bucket of (arbitrary) node -- now the closest process

	std::vector<eslocal> closest(_dMesh.esize.size());
	_eBuckets.resize(_dMesh.esize.size());

	eslocal nbegin = _nDistribution[environment->MPIrank];
	eslocal nend = _nDistribution[environment->MPIrank + 1];

	for (size_t e = 0, offset = 0; e < _dMesh.esize.size(); offset += _dMesh.esize[e++]) {
		closest[e] = _dMesh.enodes[offset];
		for (eslocal n = 1; n < _dMesh.esize[e]; n++) {
			if (nbegin <= _dMesh.enodes[offset + n] && _dMesh.enodes[offset + n] < nend) {
				if (closest[e] > _dMesh.enodes[offset + n] || closest[e] < nbegin) {
					closest[e] = _dMesh.enodes[offset + n];
				}
			} else {
				if (std::abs(closest[e] - nbegin) > std::abs(_dMesh.enodes[offset + n] - nbegin)) {
					closest[e] = _dMesh.enodes[offset + n];
				}
			}
		}
	}

	std::vector<eslocal> permutation(_dMesh.esize.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return closest[i] < closest[j]; });

	std::vector<eslocal> sNodes, rNodes;
	std::vector<uint> sBuckets, rBuckets;
	std::vector<int> targets, sources;

	sNodes.reserve(permutation.size() + 2 * environment->MPIsize);

	size_t prevsize;
	auto begin = permutation.begin();
	for (int r = 0; r < environment->MPIsize; r++) {
		prevsize = sNodes.size();
		sNodes.push_back(0);
		sNodes.push_back(r);
		sNodes.push_back(environment->MPIrank);

		auto n = begin;
		for ( ; n != permutation.end() && closest[*n] < _nDistribution[r + 1]; ++n) {
			sNodes.push_back(closest[*n]);
		}
		sNodes[prevsize] = 3 + n - begin;
		begin = n;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sNodes, rNodes)) {
		ESINFO(ERROR) << "ESPRESO internal error: ask neighbors for nodes buckets.";
	}

	std::vector<eslocal> boundaries(environment->MPIsize);
	size_t offset = 0;
	for (int r = 1; r < environment->MPIsize; r++, offset += rNodes[offset]) {
		boundaries[r] = rNodes[offset] + boundaries[r - 1];
	}
	std::sort(boundaries.begin(), boundaries.end(), [&] (eslocal i, eslocal j) {
		return rNodes[i + 2] < rNodes[j + 2];
	});

	sBuckets.reserve(rNodes.size());

	for (int r = 0; r < environment->MPIsize; r++) {
		offset = boundaries[r];
		eslocal size = rNodes[offset++] - 3;
		offset++; //skip rank
		offset++; //skip target

		sBuckets.push_back(3 + size);
		sBuckets.push_back(r);
		sBuckets.push_back(environment->MPIrank);
		for (size_t n = 0; n < size; ++n, ++offset) {
			sBuckets.push_back(_nBuckets[rNodes[offset] - _nDistribution[environment->MPIrank]]);
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuckets, rBuckets)) {
		ESINFO(ERROR) << "ESPRESO internal error: return nodes buckets.";
	}

	boundaries[0] = offset = 0;
	for (int r = 1; r < environment->MPIsize; r++, offset += rBuckets[offset]) {
		boundaries[r] = rBuckets[offset] + boundaries[r - 1];
	}
	std::sort(boundaries.begin(), boundaries.end(), [&] (eslocal i, eslocal j) {
		return rBuckets[i + 2] < rBuckets[j + 2];
	});

	size_t e = 0;
	for (int r = 0; r < environment->MPIsize; r++) {
		offset = boundaries[r];
		eslocal size = rBuckets[offset++] - 3;
		offset++; //skip rank
		offset++; //skip target

		for (size_t n = 0; n < size; ++n, ++offset) {
			_eBuckets[permutation[e++]] = rBuckets[offset];
		}
	}
}

void Converter::clusterize()
{
	double PRECISION = 0.02 * std::log2(environment->MPIsize);
	while (PRECISION * (_eDistribution.back() / environment->MPIsize) < 2) {
		PRECISION *= 2;
	}
	// allowed difference to the perfect distribution
	size_t ETOLERANCE = PRECISION * _eDistribution.back() / environment->MPIsize;

	std::vector<eslocal> edist({ 0 });
	edist.reserve(_dMesh.esize.size() + 1);
	for (size_t e = 0; e < _dMesh.esize.size(); e++) {
		edist.push_back(edist.back() + _dMesh.esize[e]);
	}

	std::vector<eslocal> npermutation(_nBuckets.size()), epermutation(_eBuckets.size());
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (eslocal i, eslocal j) { return _nBuckets[i] < _nBuckets[j]; });
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (eslocal i, eslocal j) { return _eBuckets[i] < _eBuckets[j]; });

	// PREPROCESS BUCKET SIZES
	size_t DEPTH = 1;
	while (_sfc.buckets(DEPTH++) < (size_t)environment->MPIsize);
	ESINFO(DETAILS) << "SFC DEPTH: " << DEPTH;

	size_t buckets = _sfc.buckets(DEPTH);
	uint bstep = _sfc.buckets(_sfc.depth()) / buckets;

	std::vector<std::vector<eslocal> > bucketSum(DEPTH);
	for (size_t d = 0; d < DEPTH; d++) {
		bucketSum[d].resize(_sfc.buckets(d + 1) + 1);
	}

	for (auto e = _eBuckets.begin(); e != _eBuckets.end(); ++e) {
		++bucketSum.back()[*e / bstep];
	}
	Esutils::sizesToOffsets(bucketSum.back());

	for (size_t d = DEPTH - 2; d < DEPTH; --d) {
		for (size_t b = 0; b < bucketSum[d].size(); ++b) {
			bucketSum[d][b] = bucketSum[d + 1][_sfc.bucketSize() * b];
		}
	}

	std::vector<eslocal> scounts, rcounts;
	_bucketsBorders.resize(environment->MPIsize + 1, _sfc.buckets(_sfc.depth()));
	_bucketsBorders.front() = 0;

	size_t LEVEL = 0;
	size_t bsize = _sfc.bucketSize();
	size_t coarsenig = _sfc.buckets(_sfc.depth());
	do {
		coarsenig /= bsize;
		scounts.resize(_sfc.sfcRefined(LEVEL).size() * (bsize + 1));
		rcounts.resize(_sfc.sfcRefined(LEVEL).size() * (bsize + 1));
		for (size_t b = 0, index = 0; b < _sfc.sfcRefined(LEVEL).size(); b++) {
			for (size_t i = 0; i <= bsize; i++, index++) {
				scounts[index] = bucketSum[LEVEL][bsize * _sfc.sfcRefined(LEVEL)[b] + i];
			}
		}

		MPI_Allreduce(scounts.data(), rcounts.data(), sizeof(eslocal) * scounts.size(), MPI_BYTE, MPITools::operations().sum, environment->MPICommunicator);

		_sfc.setLevel(LEVEL + 1);

		for (size_t b = 0; b < _sfc.sfcRefined(LEVEL).size(); b++) {
			size_t boffset = b * (bsize + 1);
			size_t rbegin = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), rcounts[boffset]) - _eDistribution.begin();
			size_t rend = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), rcounts[boffset + bsize]) - _eDistribution.begin();

			for (size_t r = rbegin, i = 0; r < rend; r++) {
				while (i <= bsize && rcounts[boffset + i] < _eDistribution[r] && _eDistribution[r] - rcounts[boffset + i] >= ETOLERANCE) {
					++i;
				}
				_bucketsBorders[r] = coarsenig * (bsize * _sfc.sfcRefined(LEVEL)[b] + i);
				if (rcounts[boffset + i] > _eDistribution[r] && rcounts[boffset + i] - _eDistribution[r] > ETOLERANCE) {
					_sfc.recurce(bsize * _sfc.sfcRefined(LEVEL)[b] + i - 1);
				}
			}
		}
		_sfc.finishLevel(LEVEL + 1);

	} while (++LEVEL < DEPTH && _sfc.hasLevel(LEVEL));

	// Go deeper if needed
	scounts.resize(_sfc.sfcRefined(LEVEL).size() * (bsize + 1));
	std::fill(scounts.begin(), scounts.end(), 0);
	for (size_t b = 0, index = 0; b < _sfc.sfcRefined(LEVEL).size(); b++, index += bsize + 1) {
		scounts[index] = bucketSum[LEVEL - 1][_sfc.sfcRefined(LEVEL)[b]];
	}
	rcounts.resize(scounts.size());
	std::vector<eslocal> refinedindices;

	while (LEVEL < SFCDEPTH && _sfc.hasLevel(LEVEL)) {
		coarsenig /= bsize;
		bstep /= buckets;
		for (size_t b = 0, index = 0; b < _sfc.sfcRefined(LEVEL).size(); b++, index++) {
			auto e = std::lower_bound(epermutation.begin(), epermutation.end(), coarsenig * bsize * _sfc.sfcRefined(LEVEL)[b], [&] (eslocal i, eslocal bound) { return _eBuckets[i] < bound; });
			for (size_t i = 0; i < bsize; i++, index++) {
				while (e != epermutation.end() && _eBuckets[*e] < coarsenig * bsize * _sfc.sfcRefined(LEVEL)[b] + (i + 1) * bstep) {
					++scounts[index + 1];
					++e;
				}
				scounts[index + 1] += scounts[index];
			}
		}

		MPI_Allreduce(scounts.data(), rcounts.data(), sizeof(eslocal) * scounts.size(), MPI_BYTE, MPITools::operations().sum, environment->MPICommunicator);

		_sfc.setLevel(LEVEL + 1);

		for (size_t b = 0; b < _sfc.sfcRefined(LEVEL).size(); b++) {
			size_t boffset = b * (bsize + 1);
			size_t rbegin = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), rcounts[boffset]) - _eDistribution.begin();
			size_t rend = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), rcounts[boffset + bsize]) - _eDistribution.begin();

			for (size_t r = rbegin, i = 0; r < rend; r++) {
				while (i <= bsize && rcounts[boffset + i] < _eDistribution[r] && _eDistribution[r] - rcounts[boffset + i] >= ETOLERANCE) {
					++i;
				}
				_bucketsBorders[r] = coarsenig * (bsize * _sfc.sfcRefined(LEVEL)[b] + i);
				if (rcounts[boffset + i] > _eDistribution[r] && rcounts[boffset + i] - _eDistribution[r] > ETOLERANCE) {
					_sfc.recurce(bsize * _sfc.sfcRefined(LEVEL)[b] + i - 1);
					refinedindices.push_back(boffset + i - 1);
				}
			}
		}
		_sfc.finishLevel(++LEVEL);
		Esutils::sortAndRemoveDuplicity(refinedindices);

		rcounts.swap(scounts);
		scounts.resize(_sfc.sfcRefined(LEVEL).size() * (bsize + 1));
		std::fill(scounts.begin(), scounts.end(), 0);
		for (size_t b = 0, index = 0; b < _sfc.sfcRefined(LEVEL).size(); b++, index += bsize + 1) {
			scounts[index] = rcounts[refinedindices[b]];
		}
		rcounts.resize(scounts.size());
		refinedindices.clear();
	}

	ESINFO(DETAILS) << "RECURSION DEPTH: " << LEVEL;

	std::vector<eslocal> sBuffer, rBuffer;
	sBuffer.reserve(
			5 * environment->MPIsize +
			_dMesh.esize.size() +
			_dMesh.edata.size() * sizeof(EData) / sizeof(eslocal) +
			_dMesh.enodes.size() +
			_dMesh.nIDs.size() +
			_dMesh.coordinates.size() * sizeof(Point) / sizeof(eslocal));

	size_t prevsize;
	auto nbegin = npermutation.begin();
	auto ebegin = epermutation.begin();
	for (int r = 0; r < environment->MPIsize; r++) {
		prevsize = sBuffer.size();
		sBuffer.push_back(0); // total size
		sBuffer.push_back(r); // target
		sBuffer.push_back(0); // number of elements
		sBuffer.push_back(0); // number of elements nodes
		sBuffer.push_back(0); // number of coordinates

		auto e = ebegin;
		for ( ; e != epermutation.end() && _eBuckets[*e] < _bucketsBorders[r + 1]; ++e) {
			sBuffer.push_back(_dMesh.esize[*e]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const eslocal*>(_dMesh.edata.data() + *e), reinterpret_cast<const eslocal*>(_dMesh.edata.data() + *e + 1));
			sBuffer.insert(sBuffer.end(), _dMesh.enodes.begin() + edist[*e], _dMesh.enodes.begin() + edist[*e + 1]);
			sBuffer[prevsize + 3] += edist[*e + 1] - edist[*e];
		}
		sBuffer[prevsize + 2] = e - ebegin;
		ebegin = e;

		auto n = nbegin;
		for ( ; n != npermutation.end() && _nBuckets[*n] < _bucketsBorders[r + 1]; ++n) {
			sBuffer.push_back(_dMesh.nIDs[*n]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const eslocal*>(_dMesh.coordinates.data() + *n), reinterpret_cast<const eslocal*>(_dMesh.coordinates.data() + *n + 1));
		}
		sBuffer[prevsize + 4] = n - nbegin;
		nbegin = n;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute elements according to SFC.";
	}

	_dMesh.esize.clear();
	_dMesh.edata.clear();
	_dMesh.enodes.clear();

	_dMesh.nIDs.swap(_nIDs); // keep for later usage in linkup phase
	_dMesh.coordinates.clear();

	size_t offset = 0;
	EData edata;
	Point point;
	for (int r = 0; r < environment->MPIsize; r++) {
		++offset;
		size_t esize = rBuffer[++offset];
		size_t enodes = rBuffer[++offset];
		size_t csize = rBuffer[++offset]; // coordinates
		++offset;

		for (size_t e = 0; e < esize; ++e) {
			_dMesh.esize.push_back(rBuffer[offset++]);
			memcpy(&edata, rBuffer.data() + offset, sizeof(EData));
			_dMesh.edata.push_back(edata);
			offset += sizeof(EData) / sizeof(eslocal);
			_dMesh.enodes.insert(_dMesh.enodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _dMesh.esize.back());
			offset += _dMesh.esize.back();
		}
		for (size_t c = 0; c < csize; ++c) {
			_dMesh.nIDs.push_back(rBuffer[offset]);
			++offset;
			memcpy(&point, rBuffer.data() + offset, sizeof(Point));
			_dMesh.coordinates.push_back(point);
			offset += sizeof(Point) / sizeof(eslocal);
		}
	}

	npermutation.resize(_dMesh.nIDs.size());
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.nIDs[i] < _dMesh.nIDs[j]; });

	std::vector<eslocal> newIDs;
	newIDs.reserve(npermutation.size());
	for (auto i = npermutation.begin(); i != npermutation.end(); ++i) {
		newIDs.push_back(_dMesh.nIDs[*i]);
	}
	_dMesh.nIDs.swap(newIDs);

	std::vector<Point> newCoordinates;
	newIDs.reserve(npermutation.size());
	for (auto i = npermutation.begin(); i != npermutation.end(); ++i) {
		newCoordinates.push_back(_dMesh.coordinates[*i]);
	}
	_dMesh.coordinates.swap(newCoordinates);

	_eDistribution = Communication::getDistribution(_dMesh.esize.size(), MPITools::operations().sizeToOffsetsSize_t);
}

void Converter::printSFC() {
	if (environment->MPIrank) {
		return;
	}

	std::ofstream os(Esutils::createDirectory({ Logging::outputRoot(), "VTKLEGACY_DEBUG_OUTPUT" }) + "SFC.vtk");
	os << "# vtk DataFile Version 2.0\n";
	os << "EXAMPLE\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n\n";

	size_t maxdepth = 0, n;
	while (_sfc.hasLevel(maxdepth)) {
		++maxdepth;
	}
	n = 1 << maxdepth;

	os << "POINTS " << pow(n + 1, _mesh.dimension) << " float\n";

	for (size_t k = 0; k <= n; k++) {
		for (size_t j = 0; j <= n; j++) {
			for (size_t i = 0; i <= n; i++) {
				if (_mesh.dimension == 3) {
					os << _sfc.origin().x + i * _sfc.size().x / n << " " << _sfc.origin().y + j * _sfc.size().y / n << " " << _sfc.origin().z + k * _sfc.size().z / n << " \n";
				}
				if (k == 0 && _mesh.dimension == 2) {
					os << _sfc.origin().x + i * _sfc.size().x / n << " " << _sfc.origin().y + j * _sfc.size().y / n << " 0\n";
				}
			}
		}
	}
	os << "\n";

	size_t cells = 0;
	_sfc.iterateBuckets(_bucketsBorders.front(), _bucketsBorders.back(), [&] (size_t depth, size_t index) {
		++cells;
	});

	os << "CELLS " << cells << " " << cells + _sfc.bucketSize() * cells << "\n";

	size_t bstep, x, y, z;
	if (_mesh.dimension == 2) {
		_sfc.iterateBuckets(_bucketsBorders.front(), _bucketsBorders.back(), [&] (size_t depth, size_t x, size_t y) {
			bstep = 1 << (maxdepth - depth);
			os << "4 ";
			os << (n + 1) * bstep * (y + 0) + bstep * (x + 0) << " ";
			os << (n + 1) * bstep * (y + 0) + bstep * (x + 1) << " ";
			os << (n + 1) * bstep * (y + 1) + bstep * (x + 1) << " ";
			os << (n + 1) * bstep * (y + 1) + bstep * (x + 0) << "\n";
		});
	}

	if (_mesh.dimension == 3) {
		_sfc.iterateBuckets(_bucketsBorders.front(), _bucketsBorders.back(), [&] (size_t depth, size_t x, size_t y, size_t z) {
			bstep = 1 << (maxdepth - depth);
			os << "8 ";
			os << (n + 1) * (n + 1) * bstep * (z + 0) + (n + 1) * bstep * (y + 0) + bstep * (x + 0) << " ";
			os << (n + 1) * (n + 1) * bstep * (z + 0) + (n + 1) * bstep * (y + 0) + bstep * (x + 1) << " ";
			os << (n + 1) * (n + 1) * bstep * (z + 0) + (n + 1) * bstep * (y + 1) + bstep * (x + 1) << " ";
			os << (n + 1) * (n + 1) * bstep * (z + 0) + (n + 1) * bstep * (y + 1) + bstep * (x + 0) << " ";
			os << (n + 1) * (n + 1) * bstep * (z + 1) + (n + 1) * bstep * (y + 0) + bstep * (x + 0) << " ";
			os << (n + 1) * (n + 1) * bstep * (z + 1) + (n + 1) * bstep * (y + 0) + bstep * (x + 1) << " ";
			os << (n + 1) * (n + 1) * bstep * (z + 1) + (n + 1) * bstep * (y + 1) + bstep * (x + 1) << " ";
			os << (n + 1) * (n + 1) * bstep * (z + 1) + (n + 1) * bstep * (y + 1) + bstep * (x + 0) << "\n";
		});
	}

	os << "\n";

	os << "CELL_TYPES " << cells << "\n";
	for (size_t i = 0; i < cells; i++) {
		if (_mesh.dimension == 2) {
			os << "9\n";
		}
		if (_mesh.dimension == 3) {
			os << "12\n";
		}
	}
	os << "\n";

	os << "CELL_DATA " << cells << "\n";
	os << "SCALARS MPI int 1\n";
	os << "LOOKUP_TABLE default\n";
	for (int r = 0; r < environment->MPIsize; r++) {
		_sfc.iterateBuckets(_bucketsBorders[r], _bucketsBorders[r + 1], [&] (size_t depth, size_t index) {
			os << r << "\n";
		});
	}
	os << "\n";

	os.close();
}

void Converter::linkup()
{
	// 1. Compute neighbors buckets
	// 3. Ask neighbors for coordinates
	// 4. Ask original coordinate holders for the rest nodes (for unknown nodes)
	// 5. Compute nodes neighbors
	// 6. Get real neighbors
	// 7. Re-index


	// 1. Compute neighbors buckets
	_sfc.SCFToXYZ();

	std::vector<std::pair<size_t, size_t> > neighbors;

	_sfc.iterateBuckets(_bucketsBorders[environment->MPIrank], _bucketsBorders[environment->MPIrank + 1], [&] (size_t depth, size_t index) {
		_sfc.addSFCNeighbors(depth, index, neighbors);
	});
	Esutils::sortAndRemoveDuplicity(neighbors);

	for (size_t i = 0; i < neighbors.size(); i++) {
		size_t bstep = _sfc.buckets(_sfc.depth()) / _sfc.buckets(neighbors[i].first);
		neighbors[i].first = neighbors[i].second * bstep;
		neighbors[i].second = neighbors[i].second * bstep + bstep;
	}


	if (_configuration.output.debug) {
		printSFC();
	}


	std::sort(neighbors.begin(), neighbors.end());

	size_t unique = 0;
	for (size_t i = 1; i < neighbors.size(); i++) {
		if (neighbors[i].first <= neighbors[unique].second) {
			neighbors[unique].second = neighbors[i].second;
		} else {
			neighbors[++unique] = neighbors[i];
		}
	}

	if (neighbors.size()) {
		neighbors.resize(unique + 1);
	}

	std::vector<int> nranks;
	for (size_t i = 0; i < neighbors.size(); i++) {
		int begin = std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), neighbors[i].first) - _bucketsBorders.begin();
		int end = std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), neighbors[i].second) - _bucketsBorders.begin();
		if (neighbors[i].first < _bucketsBorders[begin]) {
			--begin;
		}
		for (int r = begin; r < end; r++) {
			if (r != environment->MPIrank) {
				if (_bucketsBorders[r] != _bucketsBorders[r + 1]) {
					nranks.push_back(r);
				}
			}
		}
	}
	Esutils::sortAndRemoveDuplicity(nranks);

	// 2. Exchange elements having only one node on here

	// 3. Ask neighbors for coordinates

	// send, found, received
	std::vector<std::vector<eslocal> > sNodes(nranks.size()), fNodes(nranks.size()), rNodes(nranks.size());
	std::vector<std::vector<Point> > fCoords(nranks.size()), rCoors(nranks.size());

	std::vector<eslocal> enodes(_dMesh.enodes.begin(), _dMesh.enodes.end());
	Esutils::sortAndRemoveDuplicity(enodes);

	for (size_t id = 0, node = 0; id < _dMesh.nIDs.size() || node < enodes.size(); ++id) {
		while (node < enodes.size() && (id == _dMesh.nIDs.size() || enodes[node] < _dMesh.nIDs[id])) {
			sNodes[0].push_back(enodes[node++]);
		}
		if (node < enodes.size() && enodes[node] == _dMesh.nIDs[id]) {
			++node;
		}
	}

	for (size_t t = 1; t < nranks.size(); t++) {
		sNodes[t] = sNodes[0];
	}

	if (!Communication::exchangeUnknownSize(sNodes, rNodes, nranks)) {
		ESINFO(ERROR) << "ESPRESO internal error: request for coordinates.";
	}

	for (size_t t = 0; t < nranks.size(); t++) {
		for (size_t n = 0; n < rNodes[t].size(); n++) {
			auto node = std::lower_bound(_dMesh.nIDs.begin(), _dMesh.nIDs.end(), rNodes[t][n]);
			if (node != _dMesh.nIDs.end() && *node == rNodes[t][n]) {
				fNodes[t].push_back(*node);
				fCoords[t].push_back(_dMesh.coordinates[node - _dMesh.nIDs.begin()]);
			}
		}
	}

	if (!Communication::exchangeUnknownSize(fNodes, rNodes, nranks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested IDs.";
	}
	if (!Communication::exchangeUnknownSize(fCoords, rCoors, nranks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested coordinates.";
	}

	// 3.1 Check if all nodes are found

	size_t nodeSize = 0;
	for (size_t r = 0, i = 0; r < nranks.size(); r++) {
		nodeSize += rNodes[r].size();
	}

	// 4. If there are some unknown nodes we need to ask their original process
	// Here we use _nDistribution and _nIDs that hold old nodes distribution

	// original
	std::vector<int> oTargets, oSources;
	// send, received (that will be asked for coordinates)
	std::vector<std::vector<int> > sTargets, rTargets;
	// unknown
	std::vector<std::vector<eslocal> > uNodes;
	std::vector<std::vector<Point> > uCoords;

	if (sNodes.size() && nodeSize != sNodes.front().size()) {
		std::vector<eslocal> found, unknown(sNodes.front().size() - nodeSize);
		for (size_t r = 0, i = 0; r < nranks.size(); r++) {
			found.insert(found.end(), rNodes[r].begin(), rNodes[r].end());
		}
		Esutils::sortAndRemoveDuplicity(found);

		std::set_difference(sNodes.front().begin(), sNodes.front().end(), found.begin(), found.end(), unknown.begin());
		sNodes.clear();

		for (size_t i = 0; i < unknown.size(); i++) {
			int trank = std::lower_bound(_nDistribution.begin(), _nDistribution.end(), unknown[i] + 1) - _nDistribution.begin() - 1;
			if (oTargets.size() == 0 || oTargets.back() != trank) {
				oTargets.push_back(trank);
				sNodes.push_back({});
			}
			sNodes.back().push_back(unknown[i]);
		}
		uNodes.resize(sNodes.size());
	}

	if (!Communication::sendVariousTargets(sNodes, uNodes, oTargets, oSources)) {
		ESINFO(ERROR) << "ESPRESO internal error: request for unknown nodes.";
	}

	sTargets.resize(oSources.size());
	for (size_t t = 0; t < oSources.size(); t++) {
		for (size_t n = 0; n < uNodes[t].size(); n++) {
			auto node = std::lower_bound(_nIDs.begin(), _nIDs.end(), uNodes[t][n]);
			if (node != _dMesh.nIDs.end() && *node == uNodes[t][n]) {
				sTargets[t].push_back(std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), _nBuckets[node - _nIDs.begin()]) - _bucketsBorders.begin() - 1);
			}
		}
	}

	if (!Communication::sendVariousTargets(sTargets, rTargets, oSources)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested unknown node targets.";
	}

	sNodes.clear();

	for (size_t i = 1; i < oTargets.size(); i++) {
		uNodes[0].insert(uNodes[0].end(), uNodes[i].begin(), uNodes[i].end());
		rTargets[0].insert(rTargets[0].end(), rTargets[i].begin(), rTargets[i].end());
	}

	if (oTargets.size()) {
		oTargets.clear();
		std::vector<eslocal> upermutation(uNodes.front().size());
		std::iota(upermutation.begin(), upermutation.end(), 0);
		std::sort(upermutation.begin(), upermutation.end(), [&] (eslocal i, eslocal j) { return rTargets[0][i] < rTargets[0][i]; });

		for (size_t i = 0; i < upermutation.size(); i++) {
			if (i == 0 || rTargets[0][upermutation[i]] != rTargets[0][upermutation[i - 1]]) {
				oTargets.push_back(rTargets[0][upermutation[i]]);
				sNodes.push_back({});
			}
			sNodes.back().push_back(uNodes[0][upermutation[i]]);
		}
	}

	uNodes.clear();
	oSources.clear();

	if (!Communication::sendVariousTargets(sNodes, uNodes, oTargets, oSources)) {
		ESINFO(ERROR) << "ESPRESO internal error: request for unknown nodes.";
	}

	fCoords.clear();
	fCoords.resize(oSources.size());
	for (size_t t = 0; t < oSources.size(); t++) {
		for (size_t n = 0; n < uNodes[t].size(); n++) {
			auto node = std::lower_bound(_dMesh.nIDs.begin(), _dMesh.nIDs.end(), uNodes[t][n]);
			if (node != _dMesh.nIDs.end() && *node == uNodes[t][n]) {
				fCoords[t].push_back(_dMesh.coordinates[node - _dMesh.nIDs.begin()]);
			}
		}
	}

	if (!Communication::sendVariousTargets(fCoords, uCoords, oSources)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested unknown coordinates.";
	}

	// insert new neighbors to neighbors computed from SFC
	for (size_t i = 0; i < oTargets.size(); i++) {
		auto it = std::lower_bound(nranks.begin(), nranks.end(), oTargets[i]);
		nranks.insert(it, oTargets[i]);
		rNodes.insert(rNodes.begin() + (it - nranks.begin()), sNodes[i]);
		rCoors.insert(rCoors.begin() + (it - nranks.begin()), uCoords[i]);
		fNodes.insert(fNodes.begin() + (it - nranks.begin()), std::vector<eslocal>());
	}

	for (size_t i = 0; i < oSources.size(); i++) {
		auto it = std::lower_bound(nranks.begin(), nranks.end(), oSources[i]);
		nranks.insert(it, oSources[i]);
		rNodes.insert(rNodes.begin() + (it - nranks.begin()), std::vector<eslocal>());
		rCoors.insert(rCoors.begin() + (it - nranks.begin()), std::vector<Point>());
		fNodes.insert(fNodes.begin() + (it - nranks.begin()), rNodes[i]);
	}

	nranks.push_back(environment->MPIrank);
	std::sort(nranks.begin(), nranks.end());

	std::vector<std::vector<eslocal> > sRanks(nranks.size()), rRanks(nranks.size());


	// 5. Compute nodes neighbors
	size_t rankindex;
	std::vector<std::vector<eslocal> > nodeRequests(nranks.size());
	for (size_t r = 0, i = 0; r < nranks.size(); r++) {
		if (nranks[r] == environment->MPIrank) {
			rankindex = r;
			nodeRequests[r].swap(enodes);
		} else {
			nodeRequests[r].swap(fNodes[i++]);
		}
	}
	std::vector<eslocal> ranks, ranksOffset;
	std::vector<std::vector<eslocal>::const_iterator> rPointer(nodeRequests.size());
	for (size_t r = 0; r < nodeRequests.size(); r++) {
		rPointer[r] = std::lower_bound(nodeRequests[r].begin(), nodeRequests[r].end(), _dMesh.nIDs.front());
	}
	for (size_t n = 0; n < _dMesh.nIDs.size(); ++n) {
		ranks.clear();
		ranksOffset.clear();
		for (size_t r = 0; r < nodeRequests.size(); r++) {
			while (rPointer[r] != nodeRequests[r].end() && *rPointer[r] < _dMesh.nIDs[n]) {
				++rPointer[r];
			}
			if (rPointer[r] != nodeRequests[r].end() && *rPointer[r] == _dMesh.nIDs[n]) {
				ranksOffset.push_back(r);
				ranks.push_back(nranks[r]);
				++rPointer[r];
			}
		}
		for (size_t r = 0; r < ranks.size(); r++) {
			sRanks[ranksOffset[r]].push_back(ranksOffset.size());
			sRanks[ranksOffset[r]].insert(sRanks[ranksOffset[r]].end(), ranks.begin(), ranks.end());
		}
	}

	nodeRequests[rankindex].swap(enodes);

	// remove nodes without elements
	unique = 0;
	for (size_t id = 0, node = 0; id < _dMesh.nIDs.size(); ++id) {
		while (node < enodes.size() && enodes[node] < _dMesh.nIDs[id]) {
			++node;
		}
		if (node == enodes.size()) {
			break;
		}
		if (_dMesh.nIDs[id] == enodes[node]) {
			_dMesh.nIDs[unique] = _dMesh.nIDs[id];
			_dMesh.coordinates[unique] = _dMesh.coordinates[id];
			++unique;
		}
	}

	_dMesh.nIDs.resize(unique);
	_dMesh.coordinates.resize(unique);

	if (!Communication::exchangeUnknownSize(sRanks, rRanks, nranks)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange ranks data.";
	}

	for (size_t t = 0, i = 0; t < nranks.size(); t++) {
		if (nranks[t] != environment->MPIrank) {
			_dMesh.nIDs.insert(_dMesh.nIDs.end(), rNodes[i].begin(), rNodes[i].end());
			_dMesh.coordinates.insert(_dMesh.coordinates.end(), rCoors[i].begin(), rCoors[i].end());
			++i;
		}
	}

	// 6. Sort and re-index

	std::vector<eslocal> permutation(_dMesh.nIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _dMesh.nIDs[i] < _dMesh.nIDs[j]; });

	std::vector<std::vector<eslocal> > rankData(1), rankDistribution(1);

	rankDistribution.front().push_back(0);
	for (size_t n = 0; n < rRanks[rankindex].size(); n += rRanks[rankindex][n] + 1) {
		rankData[0].insert(rankData[0].end(), rRanks[rankindex].begin() + n + 1, rRanks[rankindex].begin() + n + 1 + rRanks[rankindex][n]);
		rankDistribution[0].push_back(rankData[0].size());
	}
	for (size_t r = 0; r < rRanks.size(); r++) {
		if (nranks[r] != environment->MPIrank) {
			for (size_t n = 0; n < rRanks[r].size(); n += rRanks[r][n] + 1) {
				rankData[0].insert(rankData[0].end(), rRanks[r].begin() + n + 1, rRanks[r].begin() + n + 1 + rRanks[r][n]);
				rankDistribution[0].push_back(rankData[0].size());
			}
		}
	}

	_mesh.nodes->size = _dMesh.nIDs.size();
	_mesh.nodes->distribution = { 0, _mesh.nodes->size };
	_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, _dMesh.nIDs);
	_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, _dMesh.coordinates);
	_mesh.nodes->ranks = new serializededata<eslocal, int>(rankDistribution, rankData);

	_mesh.nodes->permute(permutation);

	_mesh.neighbours.clear();
	Esutils::sortAndRemoveDuplicity(rankData);
	for (size_t r = 0; r < rankData[0].size(); r++) {
		if (rankData[0][r] != environment->MPIrank) {
			_mesh.neighbours.push_back(rankData[0][r]);
		}
	}
	_mesh.neighboursWithMe = _mesh.neighbours;
	_mesh.neighboursWithMe.push_back(environment->MPIrank);
	std::sort(_mesh.neighboursWithMe.begin(), _mesh.neighboursWithMe.end());

	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, _dMesh.nIDs);
}

void Converter::fillElements()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > tedist(threads);
	std::vector<std::vector<eslocal> > tnodes(threads);
	std::vector<std::vector<eslocal> > eIDs(threads), rData(threads);
	std::vector<std::vector<int> > eMat(threads), eBody(threads);
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

	_mesh.elements->dimension = _mesh.dimension;
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

	for (auto n = _mesh.elements->nodes->begin()->begin(); n != _mesh.elements->nodes->end()->begin(); ++n) {
		*n = std::lower_bound(_mesh.nodes->IDs->datatarray().begin(), _mesh.nodes->IDs->datatarray().end(), *n) - _mesh.nodes->IDs->datatarray().begin();
	}
}










