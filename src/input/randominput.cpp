
#include "randominput.h"

#include "../basis/containers/serializededata.h"
#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"
#include "../basis/logging/timeeval.h"

#include "../config/ecf/root.h"

#include "../mesh/mesh.h"
#include "../mesh/preprocessing/meshpreprocessing.h"
#include "../mesh/store/nodestore.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/boundaryregionstore.h"
#include "../old/input/loader.h"

#include "../output/result/visualization/distributedvtklegacy.h"

#include <numeric>
#include <algorithm>
#include <fstream>

using namespace espreso;

void RandomInput::buildMesh(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
{
	RandomInput(configuration, meshData, mesh);
}

RandomInput::RandomInput(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh)
: Input(configuration, meshData, mesh), _sfc(_mesh.dimension, SFCDEPTH, _meshData.coordinates)
{
	ESINFO(OVERVIEW) << "Build mesh from randomly distributed elements.";
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
	_mesh.elements->reindex(_mesh.nodes->IDs);
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Balanced loader:: elements filled.";

//	TimeEvent tpolish("polish decomposition"); tpolish.start();
//	polish();
//	tpolish.end(); timing.addEvent(tpolish);
//	ESINFO(PROGRESS2) << "Balanced loader:: decomposition polished.";

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void RandomInput::assignNBuckets()
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<size_t> cdistribution = tarray<eslocal>::distribute(threads, _meshData.nIDs.size());

	_nBuckets.resize(_meshData.coordinates.size());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t n = cdistribution[t]; n < cdistribution[t + 1]; ++n) {
			_nBuckets[n] = _sfc.getBucket(_meshData.coordinates[n]);
		}
	}
}

void RandomInput::assignEBuckets()
{ // we needs to ask a neighbor process to get bucket of (arbitrary) node -- now the closest process is chosen

	TimeEval timing("ASSIGN ELEMENTS BUCKETS");
	timing.totalTime.startWithBarrier();

	TimeEvent e1("AEB GET CLOSEST PROCESS");
	e1.start();

	std::vector<eslocal> closest(_meshData.esize.size());
	_eBuckets.resize(_meshData.esize.size());

	eslocal nbegin = _nDistribution[environment->MPIrank];
	eslocal nend = _nDistribution[environment->MPIrank + 1];

	for (size_t e = 0, offset = 0; e < _meshData.esize.size(); offset += _meshData.esize[e++]) {
		closest[e] = _meshData.enodes[offset];
		for (eslocal n = 1; n < _meshData.esize[e]; n++) {
			if (nbegin <= _meshData.enodes[offset + n] && _meshData.enodes[offset + n] < nend) {
				if (closest[e] > _meshData.enodes[offset + n] || closest[e] < nbegin) {
					closest[e] = _meshData.enodes[offset + n];
				}
			} else {
				if (std::abs(closest[e] - nbegin) > std::abs(_meshData.enodes[offset + n] - nbegin)) {
					closest[e] = _meshData.enodes[offset + n];
				}
			}
		}
	}

	e1.end();
	timing.addEvent(e1);

	TimeEvent e2("AEB PREPARE SBUFFER");
	e2.start();

	std::vector<eslocal> permutation(_meshData.esize.size());
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

	e2.end();
	timing.addEvent(e2);

	TimeEvent e3("AEB REQUEST NODES");
	e3.start();

	if (!Communication::allToAllWithDataSizeAndTarget(sNodes, rNodes)) {
		ESINFO(ERROR) << "ESPRESO internal error: ask neighbors for nodes buckets.";
	}

	e3.end();
	timing.addEvent(e3);

	TimeEvent e4("AEB COMPUTE BUCKETS");
	e4.start();

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

	e4.end();
	timing.addEvent(e4);

	TimeEvent e5("AEB RETURN BUCKETS");
	e5.start();

	if (!Communication::allToAllWithDataSizeAndTarget(sBuckets, rBuckets)) {
		ESINFO(ERROR) << "ESPRESO internal error: return nodes buckets.";
	}

	e5.end();
	timing.addEvent(e5);

	TimeEvent e6("AEB STORE BUCKETS");
	e6.start();

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

	e6.end();
	timing.addEvent(e6);

	timing.totalTime.end();
	timing.printStatsMPI();
}

void RandomInput::clusterize()
{
	TimeEval timing("CLUSTERIZE ELEMENTS");
	timing.totalTime.startWithBarrier();

	TimeEvent e1("CE PREPROCESS DATA");
	e1.start();

	double PRECISION = 0.001 * std::log2(environment->MPIsize);
	if (PRECISION * (_eDistribution.back() / environment->MPIsize) < 2) {
		PRECISION = 2.01 / (_eDistribution.back() / environment->MPIsize);
	}
	// allowed difference to the perfect distribution
	size_t ETOLERANCE = PRECISION * _eDistribution.back() / environment->MPIsize;

	std::vector<eslocal> edist({ 0 });
	edist.reserve(_meshData.esize.size() + 1);
	for (size_t e = 0; e < _meshData.esize.size(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	std::vector<eslocal> npermutation(_nBuckets.size()), epermutation(_eBuckets.size());
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (eslocal i, eslocal j) { return _nBuckets[i] < _nBuckets[j]; });
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (eslocal i, eslocal j) { return _eBuckets[i] < _eBuckets[j]; });

	e1.end();
	timing.addEvent(e1);

	// PREPROCESS BUCKET SIZES
	size_t DEPTH = 1;
	while (_sfc.buckets(DEPTH++) < (size_t)environment->MPIsize);
	ESINFO(DETAILS) << "SFC DEPTH: " << DEPTH;

	TimeEvent e2("CE COMPUTE LOCAL HISTOGRAM");
	e2.start();

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

	e2.end();
	timing.addEvent(e2);

	TimeEvent e3("CE PREPROCESSED RECURSION");
	e3.start();

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

	e3.end();
	timing.addEvent(e3);

	TimeEvent e4("CE ADHOC RECURSION");
	e4.start();

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
				while (e != epermutation.end() && _eBuckets[*e] < coarsenig * bsize * _sfc.sfcRefined(LEVEL)[b] + (i + 1) * coarsenig) {
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

	e4.end();
	timing.addEvent(e4);

	ESINFO(DETAILS) << "RECURSION DEPTH: " << LEVEL;

	TimeEvent e5("CE COMPUTE SBUFFER");
	e5.start();

	std::vector<eslocal> sBuffer, rBuffer;
	sBuffer.reserve(
			5 * environment->MPIsize +
			_meshData.esize.size() +
			_meshData.edata.size() * sizeof(EData) / sizeof(eslocal) +
			_meshData.enodes.size() +
			_meshData.nIDs.size() +
			_meshData.coordinates.size() * sizeof(Point) / sizeof(eslocal));

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
			sBuffer.push_back(_meshData.esize[*e]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const eslocal*>(_meshData.edata.data() + *e), reinterpret_cast<const eslocal*>(_meshData.edata.data() + *e + 1));
			sBuffer.insert(sBuffer.end(), _meshData.enodes.begin() + edist[*e], _meshData.enodes.begin() + edist[*e + 1]);
			sBuffer[prevsize + 3] += edist[*e + 1] - edist[*e];
		}
		sBuffer[prevsize + 2] = e - ebegin;
		ebegin = e;

		auto n = nbegin;
		for ( ; n != npermutation.end() && _nBuckets[*n] < _bucketsBorders[r + 1]; ++n) {
			sBuffer.push_back(_meshData.nIDs[*n]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const eslocal*>(_meshData.coordinates.data() + *n), reinterpret_cast<const eslocal*>(_meshData.coordinates.data() + *n + 1));
		}
		sBuffer[prevsize + 4] = n - nbegin;
		nbegin = n;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	e5.end();
	timing.addEvent(e5);

	TimeEvent e6("CE EXCHANGE DATA");
	e6.start();

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		ESINFO(ERROR) << "ESPRESO internal error: distribute elements according to SFC.";
	}

	e6.end();
	timing.addEvent(e6);

	TimeEvent e7("CE POST-PROCESS DATA");
	e7.start();

	_meshData.esize.clear();
	_meshData.edata.clear();
	_meshData.enodes.clear();

	_meshData.nIDs.swap(_nIDs); // keep for later usage in linkup phase
	_meshData.coordinates.clear();

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
			_meshData.esize.push_back(rBuffer[offset++]);
			memcpy(&edata, rBuffer.data() + offset, sizeof(EData));
			_meshData.edata.push_back(edata);
			offset += sizeof(EData) / sizeof(eslocal);
			_meshData.enodes.insert(_meshData.enodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _meshData.esize.back());
			offset += _meshData.esize.back();
		}
		for (size_t c = 0; c < csize; ++c) {
			_meshData.nIDs.push_back(rBuffer[offset]);
			++offset;
			memcpy(&point, rBuffer.data() + offset, sizeof(Point));
			_meshData.coordinates.push_back(point);
			offset += sizeof(Point) / sizeof(eslocal);
		}
	}

	npermutation.resize(_meshData.nIDs.size());
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (eslocal i, eslocal j) { return _meshData.nIDs[i] < _meshData.nIDs[j]; });

	std::vector<eslocal> newIDs;
	newIDs.reserve(npermutation.size());
	for (auto i = npermutation.begin(); i != npermutation.end(); ++i) {
		newIDs.push_back(_meshData.nIDs[*i]);
	}
	_meshData.nIDs.swap(newIDs);

	std::vector<Point> newCoordinates;
	newIDs.reserve(npermutation.size());
	for (auto i = npermutation.begin(); i != npermutation.end(); ++i) {
		newCoordinates.push_back(_meshData.coordinates[*i]);
	}
	_meshData.coordinates.swap(newCoordinates);

	_eDistribution = Communication::getDistribution(_meshData.esize.size(), MPITools::operations().sizeToOffsetsSize_t);

	e7.end();
	timing.addEvent(e7);

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void RandomInput::linkup()
{
	// 1. Compute neighbors buckets
	// 3. Ask neighbors for coordinates
	// 4. Ask original coordinate holders for the rest nodes (for unknown nodes)
	// 5. Compute nodes neighbors
	// 6. Get real neighbors
	// 7. Re-index

	size_t threads = environment->OMP_NUM_THREADS;

	TimeEval timing("LINK UP");
	timing.totalTime.startWithBarrier();

	TimeEvent e1("LU SFC to XYZ");
	e1.start();

	// 1. Compute neighbors buckets
	_sfc.SCFToXYZ();

	e1.end();
	timing.addEvent(e1);

	if (_configuration.output.debug && environment->MPIrank == 0) {
		VTKLegacyDebugInfo::spaceFillingCurve(_sfc, _bucketsBorders);
	}

	TimeEvent e2("LU GET NEIGHBORS");
	e2.start();

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

	e2.end();
	timing.addEvent(e2);

	TimeEvent e3("LU COMPUTE NODES FOR NEIGHBORS");
	e3.start();

	// 2. Exchange elements having only one node on here

	// 3. Ask neighbors for coordinates

	// send, found, received
	std::vector<std::vector<eslocal> > sNodes(nranks.size()), fNodes(nranks.size()), rNodes(nranks.size());
	std::vector<std::vector<Point> > fCoords(nranks.size()), rCoors(nranks.size());

	std::vector<eslocal> enodes(_meshData.enodes.begin(), _meshData.enodes.end());
	Esutils::sortAndRemoveDuplicity(enodes);

	for (size_t id = 0, node = 0; id < _meshData.nIDs.size() || node < enodes.size(); ++id) {
		while (node < enodes.size() && (id == _meshData.nIDs.size() || enodes[node] < _meshData.nIDs[id])) {
			sNodes[0].push_back(enodes[node++]);
		}
		if (node < enodes.size() && enodes[node] == _meshData.nIDs[id]) {
			++node;
		}
	}

	for (size_t t = 1; t < nranks.size(); t++) {
		sNodes[t] = sNodes[0];
	}

	e3.end();
	timing.addEvent(e3);

	TimeEvent e4("LU REQUEST NODES");
	e4.start();

	if (!Communication::exchangeUnknownSize(sNodes, rNodes, nranks)) {
		ESINFO(ERROR) << "ESPRESO internal error: request for coordinates.";
	}

	e4.end();
	timing.addEvent(e4);

	TimeEvent e5("LU GET COORDINATES");
	e5.start();

	for (size_t t = 0; t < nranks.size(); t++) {
		for (size_t n = 0; n < rNodes[t].size(); n++) {
			auto node = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), rNodes[t][n]);
			if (node != _meshData.nIDs.end() && *node == rNodes[t][n]) {
				fNodes[t].push_back(*node);
				fCoords[t].push_back(_meshData.coordinates[node - _meshData.nIDs.begin()]);
			}
		}
	}

	e5.end();
	timing.addEvent(e5);

	TimeEvent e6("LU RETURN COORDINATES");
	e6.start();

	if (!Communication::exchangeUnknownSize(fNodes, rNodes, nranks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested IDs.";
	}
	if (!Communication::exchangeUnknownSize(fCoords, rCoors, nranks)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested coordinates.";
	}

	e6.end();
	timing.addEvent(e6);

	TimeEvent e7("LU REQUEST FOR UNKNOWN NODES");
	e7.start();

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
			if (node != _meshData.nIDs.end() && *node == uNodes[t][n]) {
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
			auto node = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), uNodes[t][n]);
			if (node != _meshData.nIDs.end() && *node == uNodes[t][n]) {
				fCoords[t].push_back(_meshData.coordinates[node - _meshData.nIDs.begin()]);
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

	e7.end();
	timing.addEvent(e7);

	TimeEvent e8("LU COMPUTE NODES TO RANK MAP");
	e8.start();

	// 5. Compute nodes neighbors
	std::vector<std::vector<eslocal> > sRanks(nranks.size()), rRanks(nranks.size());

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
		rPointer[r] = std::lower_bound(nodeRequests[r].begin(), nodeRequests[r].end(), _meshData.nIDs.front());
	}
	for (size_t n = 0; n < _meshData.nIDs.size(); ++n) {
		ranks.clear();
		ranksOffset.clear();
		for (size_t r = 0; r < nodeRequests.size(); r++) {
			while (rPointer[r] != nodeRequests[r].end() && *rPointer[r] < _meshData.nIDs[n]) {
				++rPointer[r];
			}
			if (rPointer[r] != nodeRequests[r].end() && *rPointer[r] == _meshData.nIDs[n]) {
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
	for (size_t id = 0, node = 0; id < _meshData.nIDs.size(); ++id) {
		while (node < enodes.size() && enodes[node] < _meshData.nIDs[id]) {
			++node;
		}
		if (node == enodes.size()) {
			break;
		}
		if (_meshData.nIDs[id] == enodes[node]) {
			_meshData.nIDs[unique] = _meshData.nIDs[id];
			_meshData.coordinates[unique] = _meshData.coordinates[id];
			++unique;
		}
	}

	_meshData.nIDs.resize(unique);
	_meshData.coordinates.resize(unique);

	e8.end();
	timing.addEvent(e8);

	TimeEvent e9("LU EXCHANGE NODE TO RANK MAP");
	e9.start();

	if (!Communication::exchangeUnknownSize(sRanks, rRanks, nranks)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange ranks data.";
	}

	e9.end();
	timing.addEvent(e9);

	TimeEvent e10("LU POST-PROCESS DATA");
	e10.start();

	for (size_t t = 0, i = 0; t < nranks.size(); t++) {
		if (nranks[t] != environment->MPIrank) {
			_meshData.nIDs.insert(_meshData.nIDs.end(), rNodes[i].begin(), rNodes[i].end());
			_meshData.coordinates.insert(_meshData.coordinates.end(), rCoors[i].begin(), rCoors[i].end());
			++i;
		}
	}

	// 6. Sort and re-index

	std::vector<eslocal> permutation(_meshData.nIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return _meshData.nIDs[i] < _meshData.nIDs[j]; });

	std::vector<std::vector<eslocal> > rankData(threads), rankDistribution(threads);

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

	_mesh.nodes->size = _meshData.nIDs.size();
	_mesh.nodes->distribution = tarray<eslocal>::distribute(threads, _mesh.nodes->size);
	_mesh.nodes->IDs = new serializededata<eslocal, eslocal>(1, tarray<eslocal>(_mesh.nodes->distribution, _meshData.nIDs));
	_mesh.nodes->coordinates = new serializededata<eslocal, Point>(1, tarray<Point>(_mesh.nodes->distribution, _meshData.coordinates));


	serializededata<eslocal, eslocal>::balance(rankDistribution, rankData, &_mesh.nodes->distribution);
	_mesh.nodes->ranks = new serializededata<eslocal, int>(rankDistribution, rankData);

	_mesh.nodes->permute(permutation);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		Esutils::sortAndRemoveDuplicity(rankData[t]);
	}

	_mesh.neighboursWithMe.clear();
	for (size_t t = 0; t < threads; t++) {
		_mesh.neighboursWithMe.insert(_mesh.neighboursWithMe.end(), rankData[t].begin(), rankData[t].end());
	}
	Esutils::sortAndRemoveDuplicity(_mesh.neighboursWithMe);

	_mesh.neighbours.clear();
	for (size_t n = 0; n < _mesh.neighboursWithMe.size(); n++) {
		if (_mesh.neighboursWithMe[n] != environment->MPIrank) {
			_mesh.neighbours.push_back(_mesh.neighboursWithMe[n]);
		}
	}

	_mesh.boundaryRegions.push_back(new BoundaryRegionStore("ALL_NODES", _mesh._eclasses));
	_mesh.boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, _meshData.nIDs);

	e10.end();
	timing.addEvent(e10);

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void RandomInput::polish()
{
	size_t threads = environment->OMP_NUM_THREADS;

	_mesh.preprocessing->computeElementsCenters();

	std::vector<eslocal> partition(_mesh.elements->size);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		Point p;
		for (size_t e = _mesh.elements->distribution[t]; e < _mesh.elements->distribution[t + 1]; e++) {
			p.x = _mesh.elements->centers->datatarray()[_mesh.dimension * e + 0];
			p.y = _mesh.elements->centers->datatarray()[_mesh.dimension * e + 1];
			if (_mesh.dimension == 3) {
				p.z = _mesh.elements->centers->datatarray()[_mesh.dimension * e + 2];
			}
			partition[e] = std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), _sfc.getBucket(p) + 1) - _bucketsBorders.begin() - 1;
		}
	}

	_mesh.preprocessing->exchangeElements(partition);
}








