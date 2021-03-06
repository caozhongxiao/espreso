
#include "randominput.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/communication.h"

#include "config/ecf/root.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

#include "mesh/mesh.h"
#include "mesh/elements/element.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

void RandomInput::buildMesh(PlainMeshData &meshData, Mesh &mesh)
{
	RandomInput(meshData, mesh);
}

RandomInput::RandomInput(PlainMeshData &meshData, Mesh &mesh)
: Input(meshData, mesh), _sfc(_mesh.dimension, SFCDEPTH, _meshData.coordinates)
{
	if (info::mpi::size == 1) {
		eslog::globalerror("ESPRESO internal error: use the sequential input for building mesh on 1 MPI process.\n");
	}

	eslog::startln("MESIO: BUILD SCATTERED MESH", "MESIO");

	balance();
	eslog::checkpointln("MESIO: DATA BALANCED");

	assignRegions(_meshData.eregions, _meshData.eIDs, _eDistribution, _eregsize, _eregions);
	assignRegions(_meshData.nregions, _meshData.nIDs, _nDistribution, _nregsize, _nregions);
	eslog::checkpointln("MESIO: REGION ASSIGNED");

//	reindexRegions();

	assignNBuckets();
	eslog::checkpointln("MESIO: NODES BUCKETS COMPUTED");

	assignEBuckets();
	eslog::checkpointln("MESIO: ELEMENTS BUCKETS ASSIGNED");

	clusterize();
	eslog::checkpointln("MESIO: ELEMENTS CLUSTERED");

	sortElements();
	eslog::checkpointln("MESIO: ELEMENTS SORTED");

	linkup();
	fillNeighbors();
	eslog::checkpointln("MESIO: NEIGHBORS COMPUTED");

	sortNodes();
	eslog::checkpointln("MESIO: NODES SORTED");

	fillNodes();
	eslog::checkpointln("MESIO: NODES FILLED");

	fillElements();
	eslog::checkpointln("MESIO: ELEMENTS SORTED");

	reindexElementNodes();
	eslog::checkpointln("MESIO: ELEMENTS NODES REINDEXED");

	if (_mesh.nodes->elements == NULL) {
		_mesh.preprocessing->linkNodesAndElements();
	}

	exchangeBoundary();
	eslog::checkpointln("MESIO: BOUNDARY EXCHANGED");

	fillRegions(_meshData.eregions, _eregsize, _eregions);
	fillRegions(_meshData.nregions, _nregsize, _nregions);
	fillElementRegions();
	fillBoundaryRegions();
	fillNodeRegions();
	eslog::checkpointln("MESIO: REGIONS FILLED");

	reindexBoundaryNodes();
	eslog::endln("MESIO: BOUNDARY NODES REINDEXED");

//	polish();
}

void RandomInput::assignNBuckets()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<size_t> cdistribution = tarray<size_t>::distribute(threads, _meshData.nIDs.size());

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

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<size_t> edistribution = tarray<size_t>::distribute(threads, _meshData.esize.size());

	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	std::vector<esint> closest(_meshData.esize.size());
	_eBuckets.resize(_meshData.esize.size());

	esint nbegin = _nDistribution[info::mpi::rank];
	esint nend = _nDistribution[info::mpi::rank + 1];

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = edistribution[t], offset = _meshData._edist[edistribution[t]]; e < edistribution[t + 1]; offset += _meshData.esize[e++]) {
			closest[e] = _meshData.enodes[offset];
			for (esint n = 1; n < _meshData.esize[e]; n++) {
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
	}

	std::vector<esint> permutation(_meshData.esize.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return closest[i] < closest[j]; });

	std::vector<esint> sNodes, rNodes;
	std::vector<uint> sBuckets, rBuckets;
	std::vector<int> targets, sources;

	sNodes.reserve(permutation.size() + 2 * info::mpi::size);

	size_t prevsize;
	auto begin = permutation.begin();
	for (int r = 0; r < info::mpi::size; r++) {
		prevsize = sNodes.size();
		sNodes.push_back(0);
		sNodes.push_back(r);
		sNodes.push_back(info::mpi::rank);

		auto n = begin;
		for ( ; n != permutation.end() && closest[*n] < _nDistribution[r + 1]; ++n) {
			sNodes.push_back(closest[*n]);
		}
		sNodes[prevsize] = 3 + n - begin;
		begin = n;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sNodes, rNodes)) {
		eslog::error("ESPRESO internal error: ask neighbors for nodes buckets.\n");
	}

	std::vector<esint> boundaries(info::mpi::size);
	size_t offset = 0;
	for (int r = 1; r < info::mpi::size; r++, offset += rNodes[offset]) {
		boundaries[r] = rNodes[offset] + boundaries[r - 1];
	}
	std::sort(boundaries.begin(), boundaries.end(), [&] (esint i, esint j) {
		return rNodes[i + 2] < rNodes[j + 2];
	});

	sBuckets.reserve(rNodes.size());

	for (int r = 0; r < info::mpi::size; r++) {
		offset = boundaries[r];
		esint size = rNodes[offset++] - 3;
		offset++; //skip rank
		offset++; //skip target

		sBuckets.push_back(3 + size);
		sBuckets.push_back(r);
		sBuckets.push_back(info::mpi::rank);
		auto it = _meshData.nIDs.begin();
		for (esint n = 0; n < size; ++n, ++offset) {
			while (*it < rNodes[offset]) { ++it; }
			sBuckets.push_back(_nBuckets[it - _meshData.nIDs.begin()]);
		}
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuckets, rBuckets)) {
		eslog::error("ESPRESO internal error: return nodes buckets.\n");
	}

	boundaries[0] = offset = 0;
	for (int r = 1; r < info::mpi::size; r++, offset += rBuckets[offset]) {
		boundaries[r] = rBuckets[offset] + boundaries[r - 1];
	}
	std::sort(boundaries.begin(), boundaries.end(), [&] (esint i, esint j) {
		return rBuckets[i + 2] < rBuckets[j + 2];
	});

	size_t e = 0;
	for (int r = 0; r < info::mpi::size; r++) {
		offset = boundaries[r];
		esint size = rBuckets[offset++] - 3;
		offset++; //skip rank
		offset++; //skip target

		for (esint n = 0; n < size; ++n, ++offset) {
			_eBuckets[permutation[e++]] = rBuckets[offset];
		}
	}
}

void RandomInput::clusterize()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	esint esize = _meshData.eIDs.size();
	esize = Communication::exscan(esize);
	std::vector<esint> targetDistribution = tarray<esint>::distribute(info::mpi::size, esize);

	double PRECISION = 0.001 * std::log2(info::mpi::size);
	if (PRECISION * (targetDistribution.back() / info::mpi::size) < 2) {
		PRECISION = 2.01 / (targetDistribution.back() / info::mpi::size);
	}

	// allowed difference to the perfect distribution
	int ETOLERANCE = PRECISION * targetDistribution.back() / info::mpi::size;

	if (!_meshData._edist.size()) {
		_meshData._edist = { 0 };
		_meshData._edist.reserve(_meshData.esize.size() + 1);
		for (size_t e = 0; e < _meshData.esize.size(); e++) {
			_meshData._edist.push_back(_meshData._edist.back() + _meshData.esize[e]);
		}
	}

	std::vector<esint> npermutation(_nBuckets.size()), epermutation(_eBuckets.size());
	std::iota(npermutation.begin(), npermutation.end(), 0);
	std::sort(npermutation.begin(), npermutation.end(), [&] (esint i, esint j) { return _nBuckets[i] < _nBuckets[j]; });
	std::iota(epermutation.begin(), epermutation.end(), 0);
	std::sort(epermutation.begin(), epermutation.end(), [&] (esint i, esint j) { return _eBuckets[i] < _eBuckets[j]; });

	// PREPROCESS BUCKET SIZES
	size_t DEPTH = 2;
	while (_sfc.buckets(DEPTH++) < (size_t)info::mpi::size);

	size_t buckets = _sfc.buckets(DEPTH);
	uint bstep = _sfc.buckets(_sfc.depth()) / buckets;

	std::vector<std::vector<esint> > bucketSum(DEPTH);
	for (size_t d = 0; d < DEPTH; d++) {
		bucketSum[d].resize(_sfc.buckets(d + 1) + 1);
	}

	for (auto e = _eBuckets.begin(); e != _eBuckets.end(); ++e) {
		++bucketSum.back()[*e / bstep];
	}
	utils::sizesToOffsets(bucketSum.back());

	for (size_t d = DEPTH - 2; d < DEPTH; --d) {
		for (size_t b = 0; b < bucketSum[d].size(); ++b) {
			bucketSum[d][b] = bucketSum[d + 1][_sfc.bucketSize() * b];
		}
	}

	std::vector<esint> scounts, rcounts;
	_bucketsBorders.resize(info::mpi::size + 1, _sfc.buckets(_sfc.depth()));
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

		MPI_Allreduce(scounts.data(), rcounts.data(), sizeof(esint) * scounts.size(), MPI_BYTE, MPITools::esintOperations().sum, info::mpi::comm);

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

	// Go deeper if needed
	scounts.resize(_sfc.sfcRefined(LEVEL).size() * (bsize + 1));
	std::fill(scounts.begin(), scounts.end(), 0);
	for (size_t b = 0, index = 0; b < _sfc.sfcRefined(LEVEL).size(); b++, index += bsize + 1) {
		scounts[index] = bucketSum[LEVEL - 1][_sfc.sfcRefined(LEVEL)[b]];
	}
	rcounts.resize(scounts.size());
	std::vector<esint> refinedindices;

	while (LEVEL < SFCDEPTH && _sfc.hasLevel(LEVEL)) {
		coarsenig /= bsize;
		bstep /= buckets;
		for (size_t b = 0, index = 0; b < _sfc.sfcRefined(LEVEL).size(); b++, index++) {
			auto e = std::lower_bound(epermutation.begin(), epermutation.end(), coarsenig * bsize * _sfc.sfcRefined(LEVEL)[b], [&] (esint i, uint bound) { return _eBuckets[i] < bound; });
			for (size_t i = 0; i < bsize; i++, index++) {
				while (e != epermutation.end() && _eBuckets[*e] < coarsenig * bsize * _sfc.sfcRefined(LEVEL)[b] + (i + 1) * coarsenig) {
					++scounts[index + 1];
					++e;
				}
				scounts[index + 1] += scounts[index];
			}
		}

		MPI_Allreduce(scounts.data(), rcounts.data(), sizeof(esint) * scounts.size(), MPI_BYTE, MPITools::esintOperations().sum, info::mpi::comm);

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
		utils::sortAndRemoveDuplicity(refinedindices);

		rcounts.swap(scounts);
		scounts.resize(_sfc.sfcRefined(LEVEL).size() * (bsize + 1));
		std::fill(scounts.begin(), scounts.end(), 0);
		for (size_t b = 0, index = 0; b < _sfc.sfcRefined(LEVEL).size(); b++, index += bsize + 1) {
			scounts[index] = rcounts[refinedindices[b]];
		}
		rcounts.resize(scounts.size());
		refinedindices.clear();
	}

	_nregsize = _meshData.nregions.size() / (8 * sizeof(esint)) + 1;
	_eregsize = _meshData.eregions.size() / (8 * sizeof(esint)) + 1;

	_nregions.resize(_nregsize * _meshData.nIDs.size());
	_eregions.resize(_eregsize * _meshData.eIDs.size());

	size_t r = 0;
	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion, ++r) {
		esint byte = r / (8 * sizeof(esint));
		esint bit = 1 << (r % (8 * sizeof(esint)));

		std::vector<size_t> rdistribution = tarray<size_t>::distribute(threads, nregion->second.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = rdistribution[t]; i < rdistribution[t + 1]; ++i) {
				_nregions[_nregsize * nregion->second[i] + byte] |= bit;
			}
		}
	}
	r = 0;
	for (auto eregion = _meshData.eregions.begin(); eregion != _meshData.eregions.end(); ++eregion, ++r) {
		esint byte = r / (8 * sizeof(esint));
		esint bit = 1 << (r % (8 * sizeof(esint)));

		std::vector<size_t> rdistribution = tarray<size_t>::distribute(threads, eregion->second.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = rdistribution[t]; i < rdistribution[t + 1]; ++i) {
				_eregions[_eregsize * eregion->second[i] + byte] |= bit;
			}
		}
	}

	std::vector<esint> sBuffer, rBuffer;
	sBuffer.reserve(
			5 * info::mpi::size +
			// esize, eID, etype, body, material, regions
			(5 + _eregsize) * _meshData.esize.size() +
			_meshData.enodes.size() +
			(1 + _nregsize) * _meshData.nIDs.size() +
			_meshData.coordinates.size() * sizeof(Point) / sizeof(esint));

	size_t prevsize;
	auto nbegin = npermutation.begin();
	auto ebegin = epermutation.begin();
	for (int r = 0; r < info::mpi::size; r++) {
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
			sBuffer.insert(sBuffer.end(), _meshData.enodes.begin() + _meshData._edist[*e], _meshData.enodes.begin() + _meshData._edist[*e + 1]);
			sBuffer[prevsize + 3] += _meshData._edist[*e + 1] - _meshData._edist[*e];
		}
		sBuffer[prevsize + 2] = e - ebegin;
		ebegin = e;

		auto n = nbegin;
		for ( ; n != npermutation.end() && _nBuckets[*n] < _bucketsBorders[r + 1]; ++n) {
			sBuffer.push_back(_meshData.nIDs[*n]);
			sBuffer.insert(sBuffer.end(), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n), reinterpret_cast<const esint*>(_meshData.coordinates.data() + *n + 1));
			sBuffer.insert(sBuffer.end(), _nregions.begin() + *n * _nregsize, _nregions.begin() + (*n + 1) * _nregsize);
		}
		sBuffer[prevsize + 4] = n - nbegin;
		nbegin = n;

		sBuffer[prevsize] = sBuffer.size() - prevsize;
	}

	if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer)) {
		eslog::error("ESPRESO internal error: distribute elements according to SFC.\n");
	}

	_meshData.esize.clear();
	_meshData.eIDs.clear();
	_meshData.etype.clear();
	_meshData.body.clear();
	_meshData.material.clear();
	_meshData.enodes.clear();
	_meshData._edist.clear();
	_eregions.clear();

	_meshData.nIDs.swap(_nIDs); // keep for later usage in linkup phase
	_meshData.coordinates.clear();
	_nregions.clear();

	size_t offset = 0;
	Point point;
	for (int r = 0; r < info::mpi::size; r++) {
		++offset;
		size_t esize = rBuffer[++offset];
		++offset;
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
			offset += sizeof(Point) / sizeof(esint);
			_nregions.insert(_nregions.end(), rBuffer.begin() + offset, rBuffer.begin() + offset + _nregsize);
			offset += _nregsize;
		}
	}

	if (!_meshData.eIDs.size()) {
		eslog::error("ESPRESO internal error: a process without elements -- re-run with smaller number of MPI.\n");
	}

	auto back = _meshData.eIDs.back();
	MPI_Allgather(&back, sizeof(back), MPI_BYTE, _eDistribution.data() + 1, sizeof(back), MPI_BYTE, info::mpi::comm);
	for (size_t i = 1; i < _eDistribution.size(); i++) {
		++_eDistribution[i];
	}

	sortNodes();
}

void RandomInput::linkup()
{
	// 1. Compute neighbors buckets
	// 3. Ask neighbors for coordinates
	// 4. Ask original coordinate holders for the rest nodes (for unknown nodes)
	// 5. Compute nodes neighbors

	// 1. Compute neighbors buckets
	_sfc.SCFToXYZ();

//	VTKLegacyDebugInfo::spaceFillingCurve(_sfc, _bucketsBorders);

	std::vector<std::pair<size_t, size_t> > neighbors;

	_sfc.iterateBuckets(_bucketsBorders[info::mpi::rank], _bucketsBorders[info::mpi::rank + 1], [&] (size_t depth, size_t index) {
		_sfc.addSFCNeighbors(depth, index, neighbors);
	});
	utils::sortAndRemoveDuplicity(neighbors);

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
			if (r != info::mpi::rank) {
				if (_bucketsBorders[r] != _bucketsBorders[r + 1]) {
					_sfcNeighbors.push_back(r);
				}
			}
		}
	}
	utils::sortAndRemoveDuplicity(_sfcNeighbors);

	// 2. Exchange elements having only one node on here

	// 3. Ask neighbors for coordinates

	// send, found, received
	std::vector<std::vector<esint> > sNodes(_sfcNeighbors.size()), fNodes(_sfcNeighbors.size()), fRegions(_sfcNeighbors.size()), rNodes(_sfcNeighbors.size()), rRegions(_sfcNeighbors.size());
	std::vector<std::vector<Point> > fCoords(_sfcNeighbors.size()), rCoors(_sfcNeighbors.size());

	size_t enodesize = 0;
	size_t estart = _mesh.dimension == 3 ? 0 : 1;
	for (esint e = 0; e < _etypeDistribution[estart]; e++) {
		enodesize += _meshData.esize[e];
	}

	std::vector<esint> enodes(_meshData.enodes.begin(), _meshData.enodes.begin() + enodesize);
	utils::sortAndRemoveDuplicity(enodes);

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

	if (!Communication::exchangeUnknownSize(sNodes, rNodes, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: request for coordinates.\n");
	}

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

	if (!Communication::exchangeUnknownSize(fNodes, rNodes, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: return requested IDs.\n");
	}
	if (!Communication::exchangeUnknownSize(fRegions, rRegions, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: return requested node regions.\n");
	}
	if (!Communication::exchangeUnknownSize(fCoords, rCoors, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: return requested coordinates.\n");
	}

	// 3.1 Check if all nodes are found

	size_t nodeSize = 0;
	for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
		nodeSize += rNodes[r].size();
	}

	// 4. If there are some unknown nodes we need to ask their original process
	// Here we use _nDistribution and _nIDs that hold old nodes distribution

	// original
	std::vector<int> oTargets, oSources;
	// send, received (that will be asked for coordinates)
	std::vector<std::vector<int> > sTargets, rTargets;
	// unknown
	std::vector<std::vector<esint> > uNodes, uRegions;
	std::vector<std::vector<Point> > uCoords;

	if (sNodes.size() && nodeSize != sNodes.front().size()) {
		std::vector<esint> found, unknown(sNodes.front().size() - nodeSize);
		for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
			found.insert(found.end(), rNodes[r].begin(), rNodes[r].end());
		}
		utils::sortAndRemoveDuplicity(found);

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
		eslog::error("ESPRESO internal error: request for unknown nodes.\n");
	}

	sTargets.resize(oSources.size());
	for (size_t t = 0; t < oSources.size(); t++) {
		for (size_t n = 0; n < uNodes[t].size(); n++) {
			auto node = std::lower_bound(_nIDs.begin(), _nIDs.end(), uNodes[t][n]);
			if (node != _nIDs.end() && *node == uNodes[t][n]) {
				sTargets[t].push_back(std::lower_bound(_bucketsBorders.begin(), _bucketsBorders.end(), _nBuckets[node - _nIDs.begin()] + 1) - _bucketsBorders.begin() - 1);
			} else {
				eslog::error("ESPRESO internal error: something wrong happen during link-up phase (request for unknown node).\n");
			}
		}
	}

	if (!Communication::sendVariousTargets(sTargets, rTargets, oSources)) {
		eslog::error("ESPRESO internal error: return requested unknown node targets.\n");
	}

	uNodes.clear();

	for (size_t i = 1; i < oTargets.size(); i++) {
		sNodes[0].insert(sNodes[0].end(), sNodes[i].begin(), sNodes[i].end());
		rTargets[0].insert(rTargets[0].end(), rTargets[i].begin(), rTargets[i].end());
	}

	if (oTargets.size()) {
		oTargets.clear();
		std::vector<esint> upermutation(sNodes.front().size());
		std::iota(upermutation.begin(), upermutation.end(), 0);
		std::sort(upermutation.begin(), upermutation.end(), [&] (esint i, esint j) {
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
		eslog::error("ESPRESO internal error: request for unknown nodes.\n");
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
				eslog::error("ESPRESO internal error: something wrong happen during link-up phase.\n");
			}
		}
	}

	if (!Communication::sendVariousTargets(fRegions, uRegions, oSources)) {
		eslog::error("ESPRESO internal error: return requested unknown node regions.\n");
	}
	if (!Communication::sendVariousTargets(fCoords, uCoords, oSources)) {
		eslog::error("ESPRESO internal error: return requested unknown coordinates.\n");
	}

	// insert new neighbors to neighbors computed from SFC
	for (size_t i = 0; i < oTargets.size(); i++) {
		auto it = std::lower_bound(_sfcNeighbors.begin(), _sfcNeighbors.end(), oTargets[i]);
		size_t offset = it - _sfcNeighbors.begin();
		_sfcNeighbors.insert(it, oTargets[i]);
		rNodes.insert(rNodes.begin() + offset, sNodes[i]);
		rRegions.insert(rRegions.begin() + offset, uRegions[i]);
		rCoors.insert(rCoors.begin() + offset, uCoords[i]);
		fNodes.insert(fNodes.begin() + offset, std::vector<esint>());
	}

	for (size_t i = 0; i < oSources.size(); i++) {
		auto it = std::lower_bound(_sfcNeighbors.begin(), _sfcNeighbors.end(), oSources[i]);
		size_t offset = it - _sfcNeighbors.begin();
		if (it != _sfcNeighbors.end() && *it == oSources[i]) {
			fNodes[offset].swap(uNodes[i]);
		} else {
			_sfcNeighbors.insert(it, oSources[i]);
			rNodes.insert(rNodes.begin() + offset, std::vector<esint>());
			rRegions.insert(rRegions.begin() + offset, std::vector<esint>());
			rCoors.insert(rCoors.begin() + offset, std::vector<Point>());
			fNodes.insert(fNodes.begin() + offset, uNodes[i]);
		}
	}

	_sfcNeighbors.push_back(info::mpi::rank);
	std::sort(_sfcNeighbors.begin(), _sfcNeighbors.end());

	// 5. Compute nodes neighbors
	std::vector<std::vector<esint> > sRanks(_sfcNeighbors.size()), rRanks(_sfcNeighbors.size());

	size_t rankindex;
	std::vector<std::vector<esint> > nodeRequests(_sfcNeighbors.size());
	for (size_t r = 0, i = 0; r < _sfcNeighbors.size(); r++) {
		if (_sfcNeighbors[r] == info::mpi::rank) {
			rankindex = r;
			nodeRequests[r].swap(enodes);
		} else {
			nodeRequests[r].swap(fNodes[i++]);
		}
	}
	std::vector<esint> ranks, ranksOffset;
	std::vector<std::vector<esint>::const_iterator> rPointer(nodeRequests.size());
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

	if (!Communication::exchangeUnknownSize(sRanks, rRanks, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: exchange ranks data.\n");
	}

	for (size_t t = 0, i = 0; t < _sfcNeighbors.size(); t++) {
		if (_sfcNeighbors[t] != info::mpi::rank) {
			_meshData.nIDs.insert(_meshData.nIDs.end(), rNodes[i].begin(), rNodes[i].end());
			_nregions.insert(_nregions.end(), rRegions[i].begin(), rRegions[i].end());
			_meshData.coordinates.insert(_meshData.coordinates.end(), rCoors[i].begin(), rCoors[i].end());
			++i;
		}
	}

	size_t r = 0;
	for (auto nregion = _meshData.nregions.begin(); nregion != _meshData.nregions.end(); ++nregion, ++r) {
		esint byte = r / (8 * sizeof(esint));
		esint bit = 1 << (r % (8 * sizeof(esint)));

		nregion->second.clear();
		for (size_t i = 0; i < _meshData.nIDs.size(); ++i) {
			if (_nregions[_nregsize * i + byte] & bit) {
				nregion->second.push_back(i);
			}
		}
	}

	_meshData._nrankdist.push_back(0);
	for (size_t n = 0; n < rRanks[rankindex].size(); n += rRanks[rankindex][n] + 1) {
		_meshData._nranks.insert(_meshData._nranks.end(), rRanks[rankindex].begin() + n + 1, rRanks[rankindex].begin() + n + 1 + rRanks[rankindex][n]);
		_meshData._nrankdist.push_back(_meshData._nranks.size());
	}
	for (size_t r = 0; r < _sfcNeighbors.size(); r++) {
		if (_sfcNeighbors[r] != info::mpi::rank) {
			for (size_t n = 0; n < rRanks[r].size(); n += rRanks[r][n] + 1) {
				_meshData._nranks.insert(_meshData._nranks.end(), rRanks[r].begin() + n + 1, rRanks[r].begin() + n + 1 + rRanks[r][n]);
				_meshData._nrankdist.push_back(_meshData._nranks.size());
			}
		}
	}
}

void RandomInput::exchangeBoundary()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	size_t estart = _mesh.dimension == 3 ? 0 : 1;

	_meshData._edist.clear();
	std::vector<esint> edist = { 0 };
	edist.reserve(_meshData.eIDs.size() - _etypeDistribution[estart] + 1);
	for (esint e = 0; e < _etypeDistribution[estart]; e++) {
		edist.back() += _meshData.esize[e];
	}
	for (esint e = _etypeDistribution[estart]; e < _etypeDistribution.back(); e++) {
		edist.push_back(edist.back() + _meshData.esize[e]);
	}

	std::vector<size_t> distribution = tarray<size_t>::distribute(threads, _etypeDistribution.back() - _etypeDistribution[estart]);
	std::vector<esint> emembership(distribution.back(), -1);
	std::vector<std::vector<std::pair<esint, esint> > > etargets(threads);
	std::vector<std::vector<esint> > utargets(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> nlinks;
		int counter, known;
		std::pair<esint, esint> target;
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

			auto push = [&] (esint eID) {
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
						if (_eDistribution[info::mpi::rank] <= nlinks[i] && nlinks[i] < _eDistribution[info::mpi::rank + 1]) {
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
	std::vector<std::vector<esint> > unodes(_sfcNeighbors.size()), rnodes(_sfcNeighbors.size());
	for (size_t i = 0; i < utargets[0].size(); i++) {
		unodes[0].insert(unodes[0].end(), _meshData.enodes.begin() + edist[utargets[0][i]], _meshData.enodes.begin() + edist[utargets[0][i] + 1]);
	}

	utils::sortAndRemoveDuplicity(unodes[0]);

	for (size_t n = 1; n < _sfcNeighbors.size(); n++) {
		unodes[n] = unodes[0];
	}

	if (!Communication::exchangeUnknownSize(unodes, rnodes, _sfcNeighbors)) {
		eslog::error("ESPRESO internal error: request for unknown boundary nodes.\n");
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
		eslog::error("ESPRESO internal error: return ranks of unknown boundary nodes.\n");
	}

	std::vector<esint> found(2 * unodes[0].size(), -1);
	for (size_t n = 0; n < _sfcNeighbors.size(); n++) {
		for (size_t i = 0, noffset = 0; i < found.size(); i += 2, noffset += rLinks[n][noffset] + 1) {
			if (rLinks[n][noffset] && found[i] == -1) {
				found[i] = n;
				found[i + 1] = noffset;
			}
		}
	}

	std::vector<esint> uunodes, rrLinks;
	for (size_t i = 0; i < found.size(); i += 2) {
		if (found[i] == -1) {
			uunodes.push_back(unodes[0][i / 2]);
		}
	}

	if (!Communication::allGatherUnknownSize(uunodes)) {
		eslog::error("ESPRESO internal error: allgather unknown nodes.\n");
	}
	utils::sortAndRemoveDuplicity(uunodes);

	for (size_t i = 0; i < uunodes.size(); i++) {
		auto node = std::lower_bound(_meshData.nIDs.begin(), _meshData.nIDs.end(), uunodes[i]);
		if (node != _meshData.nIDs.end() && *node == uunodes[i]) { // i have the node
			auto ranks = _mesh.nodes->ranks->begin() + (node - _meshData.nIDs.begin());
			if (ranks->front() == info::mpi::rank) { // i am the first rank that hold the node
				auto links = _mesh.nodes->elements->cbegin() + (node - _meshData.nIDs.begin());
				rrLinks.push_back(uunodes[i]);
				rrLinks.push_back(links->size());
				rrLinks.insert(rrLinks.end(), links->begin(), links->end());
			}
		}
	}

	if (!Communication::allGatherUnknownSize(rrLinks)) {
		eslog::error("ESPRESO internal error: allgather unknown nodes links.\n");
	}

	std::vector<esint> permutation(uunodes.size());
	for (size_t i = 0, noffset = 0; i < uunodes.size(); ++i, noffset += rrLinks[noffset + 1] + 2) {
		permutation[i] = noffset;
	}
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return rrLinks[i] < rrLinks[j]; });

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

	std::vector<esint> linkDist = { 0 }, linkData;
	for (size_t i = 0; i < unodes[0].size(); i++) {
		if (found[2 * i] != -1) {
			esint lindex = found[2 * i];
			esint loffset = found[2 * i + 1];
			esint lsize = rLinks[lindex][loffset];
			linkDist.push_back(linkDist.back() + lsize);
			linkData.insert(linkData.end(), rLinks[lindex].begin() + loffset + 1, rLinks[lindex].begin() + loffset + 1 + lsize);
		} else {
			esint loffset = found[2 * i + 1];
			esint lsize = rrLinks[loffset];
			linkDist.push_back(linkDist.back() + lsize);
			linkData.insert(linkData.end(), rrLinks.begin() + loffset + 1, rrLinks.begin() + loffset + 1 + lsize);
		}
	}

	{
		std::vector<esint> nlinks;
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
						esint rank = std::lower_bound(_eDistribution.begin(), _eDistribution.end(), nlinks[i] + 1) - _eDistribution.begin() - 1;
						etargets[0].push_back(std::make_pair(rank, utargets[0][e]));
						break;
					}
				} else {
					counter = 1;
				}
			}
		}
	}



	utils::sortAndRemoveDuplicity(etargets[0]);

	std::vector<int> sRanks;
	std::vector<std::vector<esint> > sBoundary, rBoundary;

	for (size_t e = 0; e < etargets[0].size(); ++e) {
		esint eindex = etargets[0][e].second + _etypeDistribution[estart];
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
		eslog::error("ESPRESO internal error: exchange boundary elements.\n");
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
		std::vector<esint> nlinks;
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
						if (_eDistribution[info::mpi::rank] <= nlinks[i] && nlinks[i] < _eDistribution[info::mpi::rank + 1]) {
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

	std::vector<esint> esize, etype, enodes, ereg;

	for (int i = estart; i < 2; i++) {
		size_t bindex = 0;
		for (esint e = _etypeDistribution[estart]; e < _etypeDistribution.back(); ++e, ++bindex) {
			if (static_cast<int>(_mesh.edata[_meshData.etype[e]].type) == 2 - i && emembership[bindex] != -1) {
				for (auto n = edist[bindex]; n < edist[bindex + 1]; ++n) {
					enodes.push_back(_meshData.enodes[n]);
				}
				esize.push_back(_meshData.esize[e]);
				etype.push_back(_meshData.etype[e]);
				ereg.insert(ereg.end(), _eregions.begin() + _eregsize * e, _eregions.begin() + _eregsize * (e + 1));
			}
		}

		for (size_t e = _etypeDistribution.back(); e < _etypeDistribution.back() + newsize; ++e, ++bindex) {
			if (static_cast<int>(_mesh.edata[_meshData.etype[e]].type) == 2 - i && emembership[bindex] != -1) {
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
			return static_cast<int>(_mesh.edata[e].type) >= type; }) - _meshData.etype.begin()
		);
	}

	_meshData._edist.clear();
}

void RandomInput::polish()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	_mesh.preprocessing->computeElementsCenters();

	std::vector<esint> partition(_mesh.elements->size);

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








