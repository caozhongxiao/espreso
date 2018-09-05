
#include "randominput.h"

#include "../basis/containers/serializededata.h"
#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"
#include "../basis/logging/timeeval.h"

#include "../config/ecf/root.h"

#include "../mesh/mesh.h"
#include "../mesh/elements/element.h"
#include "../mesh/preprocessing/meshpreprocessing.h"
#include "../mesh/store/nodestore.h"
#include "../mesh/store/elementstore.h"
#include "../mesh/store/boundaryregionstore.h"

#include "../output/result/visualization/distributedvtklegacy.h"

#include <numeric>
#include <algorithm>
#include <fstream>

using namespace espreso;

void RandomInput::buildMesh(PlainMeshData &meshData, Mesh &mesh)
{
	RandomInput(meshData, mesh);
}

RandomInput::RandomInput(PlainMeshData &meshData, Mesh &mesh)
: Input(meshData, mesh), _sfc(_mesh.dimension, SFCDEPTH, _meshData.coordinates)
{
	if (environment->MPIsize == 1) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: use the sequential input for building mesh on 1 MPI process.";
	}

	ESINFO(OVERVIEW) << "Build mesh from randomly distributed elements.";
	TimeEval timing("Load distributed mesh");
	timing.totalTime.startWithBarrier();

	TimeEvent tdistribution("distribute mesh across processes"); tdistribution.start();
	balance();
	tdistribution.end(); timing.addEvent(tdistribution);
	ESINFO(PROGRESS2) << "Random data loader:: data balanced.";

	TimeEvent tregions("assign regions"); tregions.start();
	assignRegions(_meshData.eregions, _meshData.eIDs, _eDistribution, _eregsize, _eregions);
	assignRegions(_meshData.nregions, _meshData.nIDs, _nDistribution, _nregsize, _nregions);
	tregions.end(); timing.addEvent(tregions);
	ESINFO(PROGRESS2) << "Random data loader:: regions assigned.";

//	TimeEvent treindexregions("reindex regions"); treindexregions.start();
//	reindexRegions();
//	treindexregions.end(); timing.addEvent(treindexregions);
//	ESINFO(PROGRESS2) << "Random data loader:: regions reindexed.";

	TimeEvent tnbuckets("assign nodes to buckets"); tnbuckets.start();
	assignNBuckets();
	tnbuckets.end(); timing.addEvent(tnbuckets);
	ESINFO(PROGRESS2) << "Random data loader:: nodes buckets assigned.";

	TimeEvent tebuckets("assign elements to buckets"); tebuckets.start();
	assignEBuckets();
	tebuckets.end(); timing.addEvent(tebuckets);
	ESINFO(PROGRESS2) << "Random data loader:: elements buckets assigned.";

	TimeEvent tclusterization("clusterization"); tclusterization.start();
	clusterize();
	tclusterization.end(); timing.addEvent(tclusterization);
	ESINFO(PROGRESS2) << "Random data loader:: elements clusterized.";

	TimeEvent tesort("sort elements"); tesort.start();
	sortElements();
	tesort.end(); timing.addEvent(tesort);
	ESINFO(PROGRESS2) << "Random data loader:: elements sorted.";

	TimeEvent tlinkup("link together"); tlinkup.start();
	linkup();
	fillNeighbors();
	tlinkup.end(); timing.addEvent(tlinkup);
	ESINFO(PROGRESS2) << "Random data loader:: neighbors linked up.";

	TimeEvent tnsort("sort nodes"); tnsort.start();
	sortNodes();
	tnsort.end(); timing.addEvent(tnsort);
	ESINFO(PROGRESS2) << "Random data loader:: nodes sorted.";

	TimeEvent tnodes("fill nodes"); tnodes.start();
	fillNodes();
	tnodes.end(); timing.addEvent(tnodes);
	ESINFO(PROGRESS2) << "Random data loader:: nodes filled.";

	TimeEvent telements("fill elements"); telements.start();
	fillElements();
	telements.end(); timing.addEvent(telements);
	ESINFO(PROGRESS2) << "Random data loader:: elements filled.";

	TimeEvent treindex("reindex elements nodes"); treindex.start();
	reindexElementNodes();
	treindex.end(); timing.addEvent(treindex);
	ESINFO(PROGRESS2) << "Random data loader:: elements nodes reindexed.";

	if (_mesh.nodes->elements == NULL) {
		_mesh.preprocessing->linkNodesAndElements();
	}

	TimeEvent texchange("exchange boundary"); texchange.start();
	exchangeBoundary();
	texchange.end(); timing.addEvent(texchange);
	ESINFO(PROGRESS2) << "Random data loader:: boundary exchanged.";

	TimeEvent tboundary("fill regions"); tboundary.start();
	fillRegions(_meshData.eregions, _eregsize, _eregions);
	fillRegions(_meshData.nregions, _nregsize, _nregions);
	fillElementRegions();
	fillBoundaryRegions();
	fillNodeRegions();
	tboundary.end(); timing.addEvent(tboundary);
	ESINFO(PROGRESS2) << "Random data loader:: regions filled.";

	TimeEvent treindexbondary("reindex boundary nodes"); treindexbondary.start();
	reindexBoundaryNodes();
	treindexbondary.end(); timing.addEvent(treindexbondary);
	ESINFO(PROGRESS2) << "Random data loader:: boundary nodes reindexed.";

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
		auto it = _meshData.nIDs.begin();
		for (size_t n = 0; n < size; ++n, ++offset) {
			while (*it < rNodes[offset]) { ++it; }
			sBuckets.push_back(_nBuckets[it - _meshData.nIDs.begin()]);
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

	size_t esize = _meshData.eIDs.size();
	esize = Communication::exscan(esize);
	std::vector<size_t> targetDistribution = Esutils::getDistribution(environment->MPIsize, esize);

	double PRECISION = 0.001 * std::log2(environment->MPIsize);
	if (PRECISION * (targetDistribution.back() / environment->MPIsize) < 2) {
		PRECISION = 2.01 / (targetDistribution.back() / environment->MPIsize);
	}
	// allowed difference to the perfect distribution
	size_t ETOLERANCE = PRECISION * targetDistribution.back() / environment->MPIsize;

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
	size_t DEPTH = 2;
	while (_sfc.buckets(DEPTH++) < (size_t)environment->MPIsize);

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

		MPI_Allreduce(scounts.data(), rcounts.data(), sizeof(eslocal) * scounts.size(), MPI_BYTE, MPITools::eslocalOperations().sum, environment->MPICommunicator);

		_sfc.setLevel(LEVEL + 1);

		for (size_t b = 0; b < _sfc.sfcRefined(LEVEL).size(); b++) {
			size_t boffset = b * (bsize + 1);
			size_t rbegin = std::lower_bound(targetDistribution.begin(), targetDistribution.end(), rcounts[boffset]) - targetDistribution.begin();
			size_t rend = std::lower_bound(targetDistribution.begin(), targetDistribution.end(), rcounts[boffset + bsize]) - targetDistribution.begin();

			for (size_t r = rbegin, i = 0; r < rend; r++) {
				while (i <= bsize && rcounts[boffset + i] < targetDistribution[r] && targetDistribution[r] - rcounts[boffset + i] >= ETOLERANCE) {
					++i;
				}
				_bucketsBorders[r] = coarsenig * (bsize * _sfc.sfcRefined(LEVEL)[b] + i);
				if (rcounts[boffset + i] > targetDistribution[r] && rcounts[boffset + i] - targetDistribution[r] > ETOLERANCE) {
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

		MPI_Allreduce(scounts.data(), rcounts.data(), sizeof(eslocal) * scounts.size(), MPI_BYTE, MPITools::eslocalOperations().sum, environment->MPICommunicator);

		_sfc.setLevel(LEVEL + 1);

		for (size_t b = 0; b < _sfc.sfcRefined(LEVEL).size(); b++) {
			size_t boffset = b * (bsize + 1);
			size_t rbegin = std::lower_bound(targetDistribution.begin(), targetDistribution.end(), rcounts[boffset]) - targetDistribution.begin();
			size_t rend = std::lower_bound(targetDistribution.begin(), targetDistribution.end(), rcounts[boffset + bsize]) - targetDistribution.begin();

			for (size_t r = rbegin, i = 0; r < rend; r++) {
				while (i <= bsize && rcounts[boffset + i] < targetDistribution[r] && targetDistribution[r] - rcounts[boffset + i] >= ETOLERANCE) {
					++i;
				}
				_bucketsBorders[r] = coarsenig * (bsize * _sfc.sfcRefined(LEVEL)[b] + i);
				if (rcounts[boffset + i] > targetDistribution[r] && rcounts[boffset + i] - targetDistribution[r] > ETOLERANCE) {
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

	TimeEvent e41("CE COMPUTE REGIONDATA");
	e41.start();

	_nregsize = _meshData.nregions.size() / (8 * sizeof(eslocal)) + 1;
	_eregsize = _meshData.eregions.size() / (8 * sizeof(eslocal)) + 1;

	_nregions.resize(_nregsize * _meshData.nIDs.size());
	_eregions.resize(_eregsize * _meshData.eIDs.size());

	size_t r = 0;
	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion, ++r) {
		eslocal byte = r / (8 * sizeof(eslocal));
		eslocal bit = 1 << (r % (8 * sizeof(eslocal)));

		for (size_t i = 0; i < nregion->second.size(); ++i) {
			_nregions[_nregsize * nregion->second[i] + byte] |= bit;
		}
	}
	r = 0;
	for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion, ++r) {
		eslocal byte = r / (8 * sizeof(eslocal));
		eslocal bit = 1 << (r % (8 * sizeof(eslocal)));

		for (size_t i = 0; i < eregion->second.size(); ++i) {
			_eregions[_eregsize * eregion->second[i] + byte] |= bit;
		}
	}

	e41.end();
	timing.addEvent(e41);

	TimeEvent e5("CE COMPUTE SBUFFER");
	e5.start();

	std::vector<eslocal> sBuffer, rBuffer;
	sBuffer.reserve(
			5 * environment->MPIsize +
			// esize, eID, etype, body, material, regions
			(5 + _eregsize) * _meshData.esize.size() +
			_meshData.enodes.size() +
			(1 + _nregsize) * _meshData.nIDs.size() +
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
			sBuffer.push_back(_meshData.eIDs[*e]);
			sBuffer.push_back(_meshData.etype[*e]);
			sBuffer.push_back(_meshData.body[*e]);
			sBuffer.push_back(_meshData.material[*e]);
			sBuffer.insert(sBuffer.end(), _eregions.begin() + *e * _eregsize, _eregions.begin() + (*e + 1) * _eregsize);
			sBuffer.insert(sBuffer.end(), _meshData.enodes.begin() + edist[*e], _meshData.enodes.begin() + edist[*e + 1]);
			sBuffer[prevsize + 3] += edist[*e + 1] - edist[*e];
		}
		sBuffer[prevsize + 2] = e - ebegin;
		ebegin = e;

		auto n = nbegin;
		for ( ; n != npermutation.end() && _nBuckets[*n] < _bucketsBorders[r + 1]; ++n) {
			sBuffer.push_back(_meshData.nIDs[*n]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const eslocal*>(_meshData.coordinates.data() + *n), reinterpret_cast<const eslocal*>(_meshData.coordinates.data() + *n + 1));
			sBuffer.insert(sBuffer.end(), _nregions.begin() + *n * _nregsize, _nregions.begin() + (*n + 1) * _nregsize);
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
	_meshData.eIDs.clear();
	_meshData.etype.clear();
	_meshData.body.clear();
	_meshData.material.clear();
	_meshData.enodes.clear();
	_eregions.clear();

	_meshData.nIDs.swap(_nIDs); // keep for later usage in linkup phase
	_meshData.coordinates.clear();
	_nregions.clear();

	size_t offset = 0;
	Point point;
	for (int r = 0; r < environment->MPIsize; r++) {
		++offset;
		size_t esize = rBuffer[++offset];
		size_t enodes = rBuffer[++offset];
		size_t csize = rBuffer[++offset]; // coordinates
		++offset;

		for (size_t e = 0; e < esize; ++e) {
			_meshData.esize.push_back(rBuffer[offset++]);
			_meshData.eIDs.push_back(rBuffer[offset++]);
			_meshData.etype.push_back(rBuffer[offset++]);
			_meshData.body.push_back(rBuffer[offset++]);
			_meshData.material.push_back(rBuffer[offset++]);
			_eregions.insert(_eregions.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _eregsize);
			offset += _eregsize;
			_meshData.enodes.insert(_meshData.enodes.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _meshData.esize.back());
			offset += _meshData.esize.back();
		}
		for (size_t c = 0; c < csize; ++c) {
			_meshData.nIDs.push_back(rBuffer[offset]);
			++offset;
			memcpy(&point, rBuffer.data() + offset, sizeof(Point));
			_meshData.coordinates.push_back(point);
			offset += sizeof(Point) / sizeof(eslocal);
			_nregions.insert(_nregions.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _nregsize);
			offset += _nregsize;
		}
	}

	if (!_meshData.eIDs.size()) {
		ESINFO(ERROR) << "ESPRESO internal error: a process without elements -- re-run with smaller number of MPI.";
	}

	size_t back = _meshData.eIDs.back();
	MPI_Allgather(&back, sizeof(size_t), MPI_BYTE, _eDistribution.data() + 1, sizeof(size_t), MPI_BYTE, environment->MPICommunicator);
	for (size_t i = 1; i < _eDistribution.size(); i++) {
		++_eDistribution[i];
	}

	sortNodes();

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

//	VTKLegacyDebugInfo::spaceFillingCurve(_sfc, _bucketsBorders);

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

	for (size_t i = 0; i < neighbors.size(); i++) {
		int begin = std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), neighbors[i].first) - _bucketsBorders.begin();
		int end = std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), neighbors[i].second) - _bucketsBorders.begin();
		if (begin && neighbors[i].first < _bucketsBorders[begin]) {
			--begin;
		}
		for (int r = begin; r < end; r++) {
			if (r != environment->MPIrank) {
				if (_bucketsBorders[r] != _bucketsBorders[r + 1]) {
					_sfcNeighbors.push_back(r);
				}
			}
		}
	}
	Esutils::sortAndRemoveDuplicity(_sfcNeighbors);

	e2.end();
	timing.addEvent(e2);

	TimeEvent e3("LU COMPUTE NODES FOR NEIGHBORS");
	e3.start();

	// 2. Exchange elements having only one node on here

	// 3. Ask neighbors for coordinates

	// send, found, received
	std::vector<std::vector<eslocal> > sNodes(_sfcNeighbors.size()), fNodes(_sfcNeighbors.size()), fRegions(_sfcNeighbors.size()), rNodes(_sfcNeighbors.size()), rRegions(_sfcNeighbors.size());
	std::vector<std::vector<Point> > fCoords(_sfcNeighbors.size()), rCoors(_sfcNeighbors.size());

	size_t enodesize = 0;
	size_t estart = _mesh.dimension == 3 ? 0 : 1;
	for (size_t e = 0; e < _etypeDistribution[estart]; e++) {
		enodesize += _meshData.esize[e];
	}

	std::vector<eslocal> enodes(_meshData.enodes.begin(), _meshData.enodes.begin() + enodesize);
	Esutils::sortAndRemoveDuplicity(enodes);

	for (size_t id = 0, node = 0; id < _meshData.nIDs.size() || node < enodes.size(); ++id) {
		while (node < enodes.size() && (id == _meshData.nIDs.size() || enodes[node] < _meshData.nIDs[id])) {
			sNodes[0].push_back(enodes[node++]);
		}
		if (node < enodes.size() && enodes[node] == _meshData.nIDs[id]) {
			++node;
		}
	}

	for (size_t t = 1; t < _sfcNeighbors.size(); t++) {
		sNodes[t] = sNodes[0];
	}

	e3.end();
	timing.addEvent(e3);

	TimeEvent e4("LU REQUEST NODES");
	e4.start();

	if (!Communication::exchangeUnknownSize(sNodes, rNodes, _sfcNeighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: request for coordinates.";
	}

	e4.end();
	timing.addEvent(e4);

	TimeEvent e5("LU GET COORDINATES");
	e5.start();

	for (size_t t = 0; t < _sfcNeighbors.size(); t++) {
		auto node = _meshData.nIDs.begin();
		for (size_t n = 0; n < rNodes[t].size(); n++) {
			while (node != _meshData.nIDs.end() && *node < rNodes[t][n]) {
				++node;
			}
			if (node != _meshData.nIDs.end() && *node == rNodes[t][n]) {
				fNodes[t].push_back(*node);
				fRegions[t].insert(fRegions[t].end(), _nregions.begin() + _nregsize * (node - _meshData.nIDs.begin()), _nregions.begin() + _nregsize * (node - _meshData.nIDs.begin() + 1));
				fCoords[t].push_back(_meshData.coordinates[node - _meshData.nIDs.begin()]);
			}
		}
	}

	e5.end();
	timing.addEvent(e5);

	TimeEvent e6("LU RETURN COORDINATES");
	e6.start();

	if (!Communication::exchangeUnknownSize(fNodes, rNodes, _sfcNeighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested IDs.";
	}
	if (!Communication::exchangeUnknownSize(fRegions, rRegions, _sfcNeighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested node regions.";
	}
	if (!Communication::exchangeUnknownSize(fCoords, rCoors, _sfcNeighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested coordinates.";
	}

	e6.end();
	timing.addEvent(e6);

	TimeEvent e7("LU REQUEST FOR UNKNOWN NODES");
	e7.start();

	// 3.1 Check if all nodes are found

	size_t nodeSize = 0;
	for (size_t r = 0, i = 0; r < _sfcNeighbors.size(); r++) {
		nodeSize += rNodes[r].size();
	}

	// 4. If there are some unknown nodes we need to ask their original process
	// Here we use _nDistribution and _nIDs that hold old nodes distribution

	// original
	std::vector<int> oTargets, oSources;
	// send, received (that will be asked for coordinates)
	std::vector<std::vector<int> > sTargets, rTargets;
	// unknown
	std::vector<std::vector<eslocal> > uNodes, uRegions;
	std::vector<std::vector<Point> > uCoords;

	if (sNodes.size() && nodeSize != sNodes.front().size()) {
		std::vector<eslocal> found, unknown(sNodes.front().size() - nodeSize);
		for (size_t r = 0, i = 0; r < _sfcNeighbors.size(); r++) {
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
	} else {
		sNodes.clear();
	}

	if (!Communication::sendVariousTargets(sNodes, uNodes, oTargets, oSources)) {
		ESINFO(ERROR) << "ESPRESO internal error: request for unknown nodes.";
	}

	sTargets.resize(oSources.size());
	for (size_t t = 0; t < oSources.size(); t++) {
		for (size_t n = 0; n < uNodes[t].size(); n++) {
			auto node = std::lower_bound(_nIDs.begin(), _nIDs.end(), uNodes[t][n]);
			if (node != _nIDs.end() && *node == uNodes[t][n]) {
				sTargets[t].push_back(std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), _nBuckets[node - _nIDs.begin()] + 1) - _bucketsBorders.begin() - 1);
			} else {
				ESINFO(ERROR) << "ESPRESO internal error: something wrong happen during link-up phase (request for unknown node).";
			}
		}
	}

	if (!Communication::sendVariousTargets(sTargets, rTargets, oSources)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested unknown node targets.";
	}

	uNodes.clear();

	for (size_t i = 1; i < oTargets.size(); i++) {
		sNodes[0].insert(sNodes[0].end(), sNodes[i].begin(), sNodes[i].end());
		rTargets[0].insert(rTargets[0].end(), rTargets[i].begin(), rTargets[i].end());
	}

	if (oTargets.size()) {
		oTargets.clear();
		std::vector<eslocal> upermutation(sNodes.front().size());
		std::iota(upermutation.begin(), upermutation.end(), 0);
		std::sort(upermutation.begin(), upermutation.end(), [&] (eslocal i, eslocal j) {
			if (rTargets[0][i] == rTargets[0][j]) {
				return sNodes[0][i] < sNodes[0][j];
			}
			return rTargets[0][i] < rTargets[0][j];
		});

		for (size_t i = 0; i < upermutation.size(); i++) {
			if (i == 0 || rTargets[0][upermutation[i]] != rTargets[0][upermutation[i - 1]]) {
				oTargets.push_back(rTargets[0][upermutation[i]]);
				uNodes.push_back({});
			}
			uNodes.back().push_back(sNodes[0][upermutation[i]]);
		}
	}

	sNodes.clear();
	oSources.clear();
	sNodes.swap(uNodes);

	if (!Communication::sendVariousTargets(sNodes, uNodes, oTargets, oSources)) {
		ESINFO(ERROR) << "ESPRESO internal error: request for unknown nodes.";
	}

	fCoords.clear();
	fCoords.resize(oSources.size());
	fRegions.clear();
	fRegions.resize(oSources.size());
	for (size_t t = 0; t < oSources.size(); t++) {
		for (size_t n = 0; n < uNodes[t].size(); n++) {
			auto node = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), uNodes[t][n]);
			if (node != _meshData.nIDs.end() && *node == uNodes[t][n]) {
				fCoords[t].push_back(_meshData.coordinates[node - _meshData.nIDs.begin()]);
				fRegions[t].insert(fRegions[t].end(), _nregions.begin() + _nregsize * (node - _meshData.nIDs.begin()), _nregions.begin() + _nregsize * (node - _meshData.nIDs.begin() + 1));
			} else {
				ESINFO(ERROR) << "ESPRESO internal error: something wrong happen during link-up phase.";
			}
		}
	}

	if (!Communication::sendVariousTargets(fRegions, uRegions, oSources)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested unknown node regions.";
	}
	if (!Communication::sendVariousTargets(fCoords, uCoords, oSources)) {
		ESINFO(ERROR) << "ESPRESO internal error: return requested unknown coordinates.";
	}

	// insert new neighbors to neighbors computed from SFC
	for (size_t i = 0; i < oTargets.size(); i++) {
		auto it = std::lower_bound(_sfcNeighbors.begin(), _sfcNeighbors.end(), oTargets[i]);
		size_t offset = it - _sfcNeighbors.begin();
		_sfcNeighbors.insert(it, oTargets[i]);
		rNodes.insert(rNodes.begin() + offset, sNodes[i]);
		rRegions.insert(rRegions.begin() + offset, uRegions[i]);
		rCoors.insert(rCoors.begin() + offset, uCoords[i]);
		fNodes.insert(fNodes.begin() + offset, std::vector<eslocal>());
	}

	for (size_t i = 0; i < oSources.size(); i++) {
		auto it = std::lower_bound(_sfcNeighbors.begin(), _sfcNeighbors.end(), oSources[i]);
		size_t offset = it - _sfcNeighbors.begin();
		if (it != _sfcNeighbors.end() && *it == oSources[i]) {
			fNodes[offset].swap(uNodes[i]);
		} else {
			_sfcNeighbors.insert(it, oSources[i]);
			rNodes.insert(rNodes.begin() + offset, std::vector<eslocal>());
			rRegions.insert(rRegions.begin() + offset, std::vector<eslocal>());
			rCoors.insert(rCoors.begin() + offset, std::vector<Point>());
			fNodes.insert(fNodes.begin() + offset, uNodes[i]);
		}
	}

	_sfcNeighbors.push_back(environment->MPIrank);
	std::sort(_sfcNeighbors.begin(), _sfcNeighbors.end());

	e7.end();
	timing.addEvent(e7);

	TimeEvent e8("LU COMPUTE NODES TO RANK MAP");
	e8.start();

	// 5. Compute nodes neighbors
	std::vector<std::vector<eslocal> > sRanks(_sfcNeighbors.size()), rRanks(_sfcNeighbors.size());

	size_t rankindex;
	std::vector<std::vector<eslocal> > nodeRequests(_sfcNeighbors.size());
	for (size_t r = 0, i = 0; r < _sfcNeighbors.size(); r++) {
		if (_sfcNeighbors[r] == environment->MPIrank) {
			rankindex = r;
			nodeRequests[r].swap(enodes);
		} else {
			nodeRequests[r].swap(fNodes[i++]);
		}
	}
	std::vector<eslocal> ranks, ranksOffset;
	std::vector<std::vector<eslocal>::const_iterator> rPointer(nodeRequests.size());
	for (size_t r = 0; r < nodeRequests.size(); r++) {
		rPointer[r] = nodeRequests[r].begin();
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
				ranks.push_back(_sfcNeighbors[r]);
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
			for (size_t i = 0; i < _nregsize; i++) {
				_nregions[_nregsize * unique + i] = _nregions[_nregsize * id + i];
			}
			++unique;
		}
	}

	_meshData.nIDs.resize(unique);
	_meshData.coordinates.resize(unique);
	_nregions.resize(_nregsize * unique);

	e8.end();
	timing.addEvent(e8);

	TimeEvent e9("LU EXCHANGE NODE TO RANK MAP");
	e9.start();

	if (!Communication::exchangeUnknownSize(sRanks, rRanks, _sfcNeighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange ranks data.";
	}

	e9.end();
	timing.addEvent(e9);

	TimeEvent e10("LU POST-PROCESS DATA");
	e10.start();

	for (size_t t = 0, i = 0; t < _sfcNeighbors.size(); t++) {
		if (_sfcNeighbors[t] != environment->MPIrank) {
			_meshData.nIDs.insert(_meshData.nIDs.end(), rNodes[i].begin(), rNodes[i].end());
			_nregions.insert(_nregions.end(), rRegions[i].begin(), rRegions[i].end());
			_meshData.coordinates.insert(_meshData.coordinates.end(), rCoors[i].begin(), rCoors[i].end());
			++i;
		}
	}

	size_t r = 0;
	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion, ++r) {
		eslocal byte = r / (8 * sizeof(eslocal));
		eslocal bit = 1 << (r % (8 * sizeof(eslocal)));

		nregion->second.clear();
		for (size_t i = 0; i < _meshData.nIDs.size(); ++i) {
			if (_nregions[_nregsize * i + byte] & bit) {
				nregion->second.push_back(i);
			}
		}
	}

	_meshData.ndist.push_back(0);
	for (size_t n = 0; n < rRanks[rankindex].size(); n += rRanks[rankindex][n] + 1) {
		_meshData.nranks.insert(_meshData.nranks.end(), rRanks[rankindex].begin() + n + 1, rRanks[rankindex].begin() + n + 1 + rRanks[rankindex][n]);
		_meshData.ndist.push_back(_meshData.nranks.size());
	}
	for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
		if (_sfcNeighbors[r] != environment->MPIrank) {
			for (size_t n = 0; n < rRanks[r].size(); n += rRanks[r][n] + 1) {
				_meshData.nranks.insert(_meshData.nranks.end(), rRanks[r].begin() + n + 1, rRanks[r].begin() + n + 1 + rRanks[r][n]);
				_meshData.ndist.push_back(_meshData.nranks.size());
			}
		}
	}

	e10.end();
	timing.addEvent(e10);

	timing.totalTime.endWithBarrier();
	timing.printStatsMPI();
}

void RandomInput::exchangeBoundary()
{
	size_t threads = environment->OMP_NUM_THREADS;

	size_t estart = _mesh.dimension == 3 ? 0 : 1;

	std::vector<eslocal> edist = { 0 };
	edist.reserve(_meshData.eIDs.size() - _etypeDistribution[estart] + 1);
	for (size_t e = 0; e < _etypeDistribution[estart]; e++) {
		edist.back() += _meshData.esize[e];
	}
	for (size_t e = _etypeDistribution[estart]; e < _etypeDistribution.back(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, _etypeDistribution.back() - _etypeDistribution[estart]);
	std::vector<eslocal> emembership(distribution.back(), -1);
	std::vector<std::vector<std::pair<eslocal, eslocal> > > etargets(threads);
	std::vector<std::vector<eslocal> > utargets(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> nlinks;
		int counter, known;
		std::pair<eslocal, eslocal> target;
		for (size_t e = distribution[t]; e < distribution[t + 1]; ++e) {
			nlinks.clear();
			known = 0;
			for (auto n = edist[e]; n < edist[e + 1]; ++n) {
				auto nit = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), _meshData.enodes[n]);
				if (nit != _meshData.nIDs.end() && *nit == _meshData.enodes[n]) {
					auto links = _mesh.nodes->elements->cbegin() + (nit - _meshData.nIDs.begin());
					nlinks.insert(nlinks.end(), links->begin(), links->end());
					++known;
				}
			}
			std::sort(nlinks.begin(), nlinks.end());

			auto push = [&] (eslocal eID) {
				target.first = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), eID + 1) - _eDistribution.begin() - 1;
				target.second = e;
				etargets[t].push_back(target);
			};

			if (known == 0) {
				utargets[t].push_back(e);
				continue;
			}

			counter = 1;
			if (known == 1) {
				push(nlinks.front());
			}
			for (size_t i = 1; i < nlinks.size(); ++i) {
				if (nlinks[i - 1] == nlinks[i]) {
					++counter;
					if (counter == edist[e + 1] - edist[e]) {
						if (_eDistribution[environment->MPIrank] <= nlinks[i] && nlinks[i] < _eDistribution[environment->MPIrank + 1]) {
							emembership[e] = nlinks[i];
						} else {
							push(nlinks[i]);
						}
						break;
					}
				} else {
					counter = 1;
				}
				if (counter == known) {
					push(nlinks[i]);
				}
			}
		}
	}

	for (size_t t = 1; t < threads; t++) {
		etargets[0].insert(etargets[0].end(), etargets[t].begin(), etargets[t].end());
		utargets[0].insert(utargets[0].end(), utargets[t].begin(), utargets[t].end());
	}


	/// FIND TARGETS FOR FACES WITH UNKNOWN NODES

	std::vector<std::vector<int> > fLinks(_sfcNeighbors.size()), rLinks(_sfcNeighbors.size());
	std::vector<std::vector<eslocal> > unodes(_sfcNeighbors.size()), rnodes(_sfcNeighbors.size());
	for (size_t i = 0; i < utargets[0].size(); i++) {
		unodes[0].insert(unodes[0].end(), _meshData.enodes.begin() + edist[utargets[0][i]], _meshData.enodes.begin() + edist[utargets[0][i] + 1]);
	}

	Esutils::sortAndRemoveDuplicity(unodes[0]);

	for (size_t n = 1; n < _sfcNeighbors.size(); n++) {
		unodes[n] = unodes[0];
	}

	if (!Communication::exchangeUnknownSize(unodes, rnodes, _sfcNeighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: request for unknown boundary nodes.";
	}

	for (size_t t = 0; t < _sfcNeighbors.size(); t++) {
		auto node = _meshData.nIDs.begin();
		for (size_t n = 0; n < rnodes[t].size(); n++) {
			while (node != _meshData.nIDs.end() && *node < rnodes[t][n]) {
				++node;
			}
			if (node != _meshData.nIDs.end() && *node == rnodes[t][n]) {
				auto links = _mesh.nodes->elements->cbegin() + (node - _meshData.nIDs.begin());
				fLinks[t].push_back(links->size());
				fLinks[t].insert(fLinks[t].end(), links->begin(), links->end());
			} else {
				fLinks[t].push_back(0);
			}

		}
	}

	if (!Communication::exchangeUnknownSize(fLinks, rLinks, _sfcNeighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: return ranks of unknown boundary nodes.";
	}

	std::vector<eslocal> found(2 * unodes[0].size(), -1);
	for (size_t n = 0; n < _sfcNeighbors.size(); n++) {
		for (size_t i = 0, noffset = 0; i < found.size(); i += 2, noffset += rLinks[n][noffset] + 1) {
			if (rLinks[n][noffset] && found[i] == -1) {
				found[i] = n;
				found[i + 1] = noffset;
			}
		}
	}

	std::vector<eslocal> uunodes, rrLinks;
	for (size_t i = 0; i < found.size(); i += 2) {
		if (found[i] == -1) {
			uunodes.push_back(unodes[0][i / 2]);
		}
	}

	if (!Communication::allGatherUnknownSize(uunodes)) {
		ESINFO(ERROR) << "ESPRESO internal error: allgather unknown nodes.";
	}
	Esutils::sortAndRemoveDuplicity(uunodes);

	for (size_t i = 0; i < uunodes.size(); i++) {
		auto node = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), uunodes[i]);
		if (node != _meshData.nIDs.end() && *node == uunodes[i]) { // i have the node
			auto ranks = _mesh.nodes->ranks->begin() + (node - _meshData.nIDs.begin());
			if (ranks->front() == environment->MPIrank) { // i am the first rank that hold the node
				auto links = _mesh.nodes->elements->cbegin() + (node - _meshData.nIDs.begin());
				rrLinks.push_back(uunodes[i]);
				rrLinks.push_back(links->size());
				rrLinks.insert(rrLinks.end(), links->begin(), links->end());
			}
		}
	}

	if (!Communication::allGatherUnknownSize(rrLinks)) {
		ESINFO(ERROR) << "ESPRESO internal error: allgather unknown nodes links.";
	}

	std::vector<eslocal> permutation(uunodes.size());
	for (size_t i = 0, noffset = 0; i < uunodes.size(); ++i, noffset += rrLinks[noffset + 1] + 2) {
		permutation[i] = noffset;
	}
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return rrLinks[i] < rrLinks[j]; });

	for (size_t i = 0, j = 0; i < unodes[0].size(); ++i) {
		while (j < uunodes.size() && uunodes[j] < unodes[0][i]) {
			++j;
		}
		if (j < uunodes.size() && unodes[0][i] == uunodes[j]) {
			if (found[2 * i] == -1) {
				found[2 * i + 1] = permutation[j] + 1;
			}
			++j;
		}
	}

	std::vector<eslocal> linkDist = { 0 }, linkData;
	for (size_t i = 0; i < unodes[0].size(); i++) {
		if (found[2 * i] != -1) {
			eslocal lindex = found[2 * i];
			eslocal loffset = found[2 * i + 1];
			eslocal lsize = rLinks[lindex][loffset];
			linkDist.push_back(linkDist.back() + lsize);
			linkData.insert(linkData.end(), rLinks[lindex].begin() + loffset + 1, rLinks[lindex].begin() + loffset + 1 + lsize);
		} else {
			eslocal loffset = found[2 * i + 1];
			eslocal lsize = rrLinks[loffset];
			linkDist.push_back(linkDist.back() + lsize);
			linkData.insert(linkData.end(), rrLinks.begin() + loffset + 1, rrLinks.begin() + loffset + 1 + lsize);
		}
	}

	{
		std::vector<eslocal> nlinks;
		int counter;
		for (size_t e = 0; e < utargets[0].size(); ++e) {
			nlinks.clear();
			for (auto n = edist[utargets[0][e]]; n < edist[utargets[0][e] + 1]; ++n) {
				auto nit = std::lower_bound(unodes[0].begin(), unodes[0].end(), _meshData.enodes[n]);
				if (nit != unodes[0].end() && *nit == _meshData.enodes[n]) {
					nlinks.insert(nlinks.end(), linkData.begin() + linkDist[nit - unodes[0].begin()], linkData.begin() + linkDist[nit - unodes[0].begin() + 1]);
				}
			}
			std::sort(nlinks.begin(), nlinks.end());

			counter = 1;
			for (size_t i = 1; i < nlinks.size(); ++i) {
				if (nlinks[i - 1] == nlinks[i]) {
					++counter;
					if (counter == edist[utargets[0][e] + 1] - edist[utargets[0][e]]) {
						eslocal rank = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), nlinks[i] + 1) - _eDistribution.begin() - 1;
						etargets[0].push_back(std::make_pair(rank, utargets[0][e]));
						break;
					}
				} else {
					counter = 1;
				}
			}
		}
	}



	Esutils::sortAndRemoveDuplicity(etargets[0]);

	std::vector<int> sRanks;
	std::vector<std::vector<eslocal> > sBoundary, rBoundary;

	for (size_t e = 0; e < etargets[0].size(); ++e) {
		eslocal eindex = etargets[0][e].second + _etypeDistribution[estart];
		if (!sRanks.size() || sRanks.back() != etargets[0][e].first) {
			sRanks.push_back(etargets[0][e].first);
			sBoundary.push_back({});
		}
		sBoundary.back().push_back(_meshData.esize[eindex]);
		sBoundary.back().push_back(_meshData.etype[eindex]);
		sBoundary.back().insert(sBoundary.back().end(), _eregions.begin() + eindex * _eregsize, _eregions.begin() + (eindex + 1) * _eregsize);
		sBoundary.back().insert(sBoundary.back().end(), _meshData.enodes.begin() + edist[etargets[0][e].second], _meshData.enodes.begin() + edist[etargets[0][e].second + 1]);
	}

	if (!Communication::sendVariousTargets(sBoundary, rBoundary, sRanks)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange boundary elements.";
	}

	for (size_t r = 1; r < rBoundary.size(); r++) {
		rBoundary[0].insert(rBoundary[0].end(), rBoundary[r].begin(), rBoundary[r].end());
	}

	size_t newsize = 0;
	for (size_t i = 0; rBoundary.size() && i < rBoundary[0].size(); ++newsize) {
		_meshData.esize.push_back(rBoundary[0][i++]);
		edist.push_back(edist.back() + _meshData.esize.back());
		_meshData.etype.push_back(rBoundary[0][i++]);
		_eregions.insert(_eregions.end(), rBoundary[0].begin() + i, rBoundary[0].begin() + i + _eregsize);
		i += _eregsize;
		_meshData.enodes.insert(_meshData.enodes.end(), rBoundary[0].begin() + i, rBoundary[0].begin() + i + _meshData.esize.back());
		i += _meshData.esize.back();
	}

	emembership.resize(emembership.size() + newsize, -1);

	{
		std::vector<eslocal> nlinks;
		int counter;
		for (size_t e = distribution.back(); e < distribution.back() + newsize; ++e) {
			nlinks.clear();
			for (auto n = edist[e]; n < edist[e + 1]; ++n) {
				auto nit = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), _meshData.enodes[n]);
				if (nit != _meshData.nIDs.end() && *nit == _meshData.enodes[n]) {
					auto links = _mesh.nodes->elements->cbegin() + (nit - _meshData.nIDs.begin());
					nlinks.insert(nlinks.end(), links->begin(), links->end());
				}
			}
			std::sort(nlinks.begin(), nlinks.end());
			counter = 1;
			for (size_t i = 1; i < nlinks.size(); ++i) {
				if (nlinks[i - 1] == nlinks[i]) {
					++counter;
					if (counter == edist[e + 1] - edist[e]) {
						if (_eDistribution[environment->MPIrank] <= nlinks[i] && nlinks[i] < _eDistribution[environment->MPIrank + 1]) {
							emembership[e] = nlinks[i];
						}
						break;
					}
				} else {
					counter = 1;
				}
			}
		}
	}

	std::vector<eslocal> esize, etype, enodes, ereg;

	for (int i = estart; i < 2; i++) {
		size_t bindex = 0;
		for (size_t e = _etypeDistribution[estart]; e < _etypeDistribution.back(); ++e, ++bindex) {
			if (static_cast<int>(_mesh._eclasses[0][_meshData.etype[e]].type) == 2 - i && emembership[bindex] != -1) {
				for (auto n = edist[bindex]; n < edist[bindex + 1]; ++n) {
					enodes.push_back(_meshData.enodes[n]);
				}
				esize.push_back(_meshData.esize[e]);
				etype.push_back(_meshData.etype[e]);
				ereg.insert(ereg.end(), _eregions.begin() + _eregsize * e, _eregions.begin() + _eregsize * (e + 1));
			}
		}

		for (size_t e = _etypeDistribution.back(); e < _etypeDistribution.back() + newsize; ++e, ++bindex) {
			if (static_cast<int>(_mesh._eclasses[0][_meshData.etype[e]].type) == 2 - i && emembership[bindex] != -1) {
				for (auto n = edist[bindex]; n < edist[bindex + 1]; ++n) {
					enodes.push_back(_meshData.enodes[n]);
				}
				esize.push_back(_meshData.esize[e]);
				etype.push_back(_meshData.etype[e]);
				ereg.insert(ereg.end(), _eregions.begin() + _eregsize * e, _eregions.begin() + _eregsize * (e + 1));
			}
		}
	}

	_meshData.esize.resize(_etypeDistribution[estart]);
	_meshData.etype.resize(_etypeDistribution[estart]);
	_meshData.enodes.resize(edist.front());
	_eregions.resize(_etypeDistribution[estart] * _eregsize);

	_meshData.esize.insert(_meshData.esize.end(), esize.begin(), esize.end());
	_meshData.etype.insert(_meshData.etype.end(), etype.begin(), etype.end());
	_meshData.enodes.insert(_meshData.enodes.end(), enodes.begin(), enodes.end());
	_eregions.insert(_eregions.end(), ereg.begin(), ereg.end());

	_etypeDistribution.clear();
	for (int type = static_cast<int>(Element::TYPE::VOLUME); type > static_cast<int>(Element::TYPE::POINT); --type) {
		_etypeDistribution.push_back(std::lower_bound(_meshData.etype.begin(), _meshData.etype.end(), type, [&] (int e, int type) {
			return static_cast<int>(_mesh._eclasses[0][e].type) >= type; }) - _meshData.etype.begin()
		);
	}
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








