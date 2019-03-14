
#include "meshpreprocessing.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

#include "mesh/mesh.h"

#include "mesh/elements/element.h"
#include "mesh/store/store.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"

#include <algorithm>
#include <numeric>

#include "wrappers/metis/metiswrapper.h"
#include "wrappers/metis/parmetiswrapper.h"

using namespace espreso;

void MeshPreprocessing::arrangeNodes()
{
	if (_mesh->elements->domainDistribution.size() == 0) {
		return;
	}

	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	std::vector<esint> externalBoundary, internalBoundary;
	this->computeBoundaryNodes(externalBoundary, internalBoundary);

	eslog::startln("MESH: ARRANGE NODES", "ARRANGE NODES");

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<esint> eDist = _mesh->elements->gatherElementsDistribution();
	std::vector<esint> dProcDist = _mesh->elements->gatherDomainsProcDistribution();
	std::vector<std::vector<esint> > domainsDistribution(threads), domainsData(threads);
	std::vector<std::vector<int> > domainsProcs(threads);

	domainsDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto elems = _mesh->nodes->elements->cbegin(t); elems != _mesh->nodes->elements->cend(t); ++elems) {
			auto prev = eDist.begin() - 1;
			auto domain = eDist.begin();
			for (auto e = elems->begin(); e != elems->end(); ++e) {
				domain = std::lower_bound(domain, eDist.end(), *e + 1) - 1;
				if (prev != domain) {
					domainsData[t].push_back(domain - eDist.begin());
					domainsProcs[t].push_back(std::lower_bound(dProcDist.begin(), dProcDist.end(), domainsData[t].back() + 1) - dProcDist.begin() - 1);
				}
				prev = domain;
			}

			domainsDistribution[t].push_back(domainsData[t].size());
		}
	}
	utils::threadDistributionToFullDistribution(domainsDistribution);

	for (size_t t = 1; t < threads; t++) {
		domainsDistribution[0].insert(domainsDistribution[0].end(), domainsDistribution[t].begin(), domainsDistribution[t].end());
		domainsData[0].insert(domainsData[0].end(), domainsData[t].begin(), domainsData[t].end());
		domainsProcs[0].insert(domainsProcs[0].end(), domainsProcs[t].begin(), domainsProcs[t].end());
	}

	std::vector<esint> permutation;
	permutation.reserve(_mesh->nodes->size);
	permutation.insert(permutation.end(), externalBoundary.begin(), externalBoundary.end());
	permutation.insert(permutation.end(), internalBoundary.begin(), internalBoundary.end());

	auto ebpointer = externalBoundary.begin();
	auto ibpointer = internalBoundary.begin();
	for (esint n = 0; n < _mesh->nodes->size; ++n) {
		if (ebpointer != externalBoundary.end() && *ebpointer == n) {
			++ebpointer;
			continue;
		}
		if (ibpointer != internalBoundary.end() && *ibpointer == n) {
			++ibpointer;
			continue;
		}
		permutation.push_back(n);
	}

	auto comp = [&] (esint i, esint j) {
		esint di = domainsDistribution[0][i], isize = domainsDistribution[0][i + 1] - domainsDistribution[0][i];
		esint dj = domainsDistribution[0][j], jsize = domainsDistribution[0][j + 1] - domainsDistribution[0][j];

		if (isize == jsize) {
			for (esint d = 0; d < isize; d++) {
				if (domainsData[0][di + d] != domainsData[0][dj + d]) {
					return domainsData[0][di + d] < domainsData[0][dj + d];
				}
			}
		}
		return isize > jsize;
	};

	std::sort(permutation.begin(), permutation.begin() + externalBoundary.size(), comp);
	std::sort(permutation.begin() + externalBoundary.size(), permutation.begin() + externalBoundary.size() + internalBoundary.size(), comp);
	std::sort(permutation.begin() + externalBoundary.size() + internalBoundary.size(), permutation.end(), comp);

	auto equalNeighs = [&] (esint i, esint j) {
		esint di = domainsDistribution[0][i], isize = domainsDistribution[0][i + 1] - domainsDistribution[0][i];
		esint dj = domainsDistribution[0][j], jsize = domainsDistribution[0][j + 1] - domainsDistribution[0][j];
		if (isize != jsize) {
			return false;
		}
		for (esint n = 0; n < isize; ++n) {
			if (domainsData[0][di + n] != domainsData[0][dj + n]) {
				return false;
			}
		}
		return true;
	};

	std::vector<std::vector<esint> > iBoundary(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = _mesh->nodes->distribution[t]; i < _mesh->nodes->distribution[t + 1]; ++i) {
			if (i > 0 && !equalNeighs(permutation[i], permutation[i - 1])) {
				iBoundary[t].push_back(i);
			}
		}
	}
	if (externalBoundary.size()) {
		iBoundary[0].push_back(externalBoundary.size());
	}
	if (internalBoundary.size()) {
		iBoundary[0].push_back(externalBoundary.size() + internalBoundary.size());
	}

	utils::mergeThreadedUniqueData(iBoundary);
	if (iBoundary[0].back() == _mesh->nodes->size) {
		iBoundary[0].pop_back();
	}

	std::vector<esint> intervalDomainsDistribution, intervalDomainsData;
	std::vector<int> intervalDomainsProcs;

	intervalDomainsDistribution.push_back(0);
	_mesh->nodes->pintervals.push_back(ProcessInterval(0, 0));
	intervalDomainsData.insert(intervalDomainsData.end(), domainsData[0].begin() + domainsDistribution[0][permutation[0]], domainsData[0].begin() + domainsDistribution[0][permutation[0] + 1]);
	intervalDomainsProcs.insert(intervalDomainsProcs.end(), domainsProcs[0].begin() + domainsDistribution[0][permutation[0]], domainsProcs[0].begin() + domainsDistribution[0][permutation[0] + 1]);
	intervalDomainsDistribution.push_back(intervalDomainsData.size());
	for (size_t i = 0; i < iBoundary[0].size(); ++i) {
		_mesh->nodes->pintervals.back().end = iBoundary[0][i];
		_mesh->nodes->pintervals.push_back(ProcessInterval(iBoundary[0][i], 0));
		intervalDomainsData.insert(intervalDomainsData.end(), domainsData[0].begin() + domainsDistribution[0][permutation[iBoundary[0][i]]], domainsData[0].begin() + domainsDistribution[0][permutation[iBoundary[0][i]] + 1]);
		intervalDomainsProcs.insert(intervalDomainsProcs.end(), domainsProcs[0].begin() + domainsDistribution[0][permutation[iBoundary[0][i]]], domainsProcs[0].begin() + domainsDistribution[0][permutation[iBoundary[0][i]] + 1]);
		intervalDomainsDistribution.push_back(intervalDomainsData.size());
	}
	_mesh->nodes->pintervals.back().end = _mesh->nodes->size;
	_mesh->nodes->idomains = new serializededata<esint, esint>(intervalDomainsDistribution, intervalDomainsData);
	_mesh->nodes->iranks = new serializededata<esint, int>(intervalDomainsDistribution, intervalDomainsProcs);

	esint externalIntervals = 0;
	auto iti = _mesh->nodes->pintervals.begin();
	while (iti != _mesh->nodes->pintervals.end() && iti->end <= (esint)externalBoundary.size()) {
		++externalIntervals;
		++iti;
	}

	#pragma omp parallel for
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i) {
		std::sort(permutation.begin() + _mesh->nodes->pintervals[i].begin, permutation.begin() + _mesh->nodes->pintervals[i].end, [&] (esint i, esint j) {
			return _mesh->nodes->IDs->datatarray()[i] < _mesh->nodes->IDs->datatarray()[j];
		});
	}

	std::vector<esint> ipermutation(_mesh->nodes->pintervals.size());
	std::iota(ipermutation.begin(), ipermutation.end(), 0);
	std::sort(ipermutation.begin(), ipermutation.end(), [&] (esint i, esint j) {
		auto n1 = _mesh->nodes->idomains->cbegin() + i;
		auto n2 = _mesh->nodes->idomains->cbegin() + j;

		esint size = std::min(n1->size(), n2->size());
		for (esint i = 0; i < size; ++i) {
			if ((*n1)[i] != (*n2)[i]) {
				return (*n1)[i] < (*n2)[i];
			}
		}
		if (n1->size() == n2->size()) {
			return i < j;
		}
		return n1->size() > n2->size();
	});

	_mesh->nodes->idomains->permute(ipermutation, _mesh->nodes->idomains->boundarytarray().distribution());
	_mesh->nodes->iranks->permute(ipermutation, _mesh->nodes->iranks->boundarytarray().distribution());
	std::vector<ProcessInterval> permutedIntervals;
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); i++) {
		permutedIntervals.push_back(_mesh->nodes->pintervals[ipermutation[i]]);
		if (ipermutation[i] < externalIntervals) {
			_mesh->nodes->externalIntervals.push_back(i);
		}
	}

	std::vector<esint> finalpermutation;
	finalpermutation.reserve(permutation.size());

	_mesh->nodes->pintervals = permutedIntervals;
	esint ioffset = 0;
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i) {
		finalpermutation.insert(finalpermutation.end(), permutation.begin() + _mesh->nodes->pintervals[i].begin, permutation.begin() + _mesh->nodes->pintervals[i].end);
		esint isize = _mesh->nodes->pintervals[i].end - _mesh->nodes->pintervals[i].begin;
		_mesh->nodes->pintervals[i].begin = ioffset;
		ioffset += isize;
		_mesh->nodes->pintervals[i].end = ioffset;
	}

	_mesh->nodes->dintervals.resize(_mesh->elements->ndomains);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (esint d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			esint doffset = 0;
			auto domains = _mesh->nodes->idomains->cbegin();
			for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i, ++domains) {
				auto iit = std::lower_bound(domains->begin(), domains->end(), _mesh->elements->firstDomain + d);
				if (iit != domains->end() && *iit == _mesh->elements->firstDomain + d) {
					_mesh->nodes->dintervals[d].push_back(DomainInterval(_mesh->nodes->pintervals[i].begin, _mesh->nodes->pintervals[i].end, i, doffset));
					doffset += _mesh->nodes->pintervals[i].end - _mesh->nodes->pintervals[i].begin;
				}
			}
		}
	}

	auto n2i = [ & ] (int neighbour) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbour) - _mesh->neighbours.begin();
	};

	std::vector<std::vector<esint> > sOffset(_mesh->neighbours.size()), rOffset(_mesh->neighbours.size());
	auto ranks = _mesh->nodes->iranks->cbegin();
	esint goffset = 0;
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i, ++ranks) {
		_mesh->nodes->pintervals[i].sourceProcess = ranks->front();
		if (ranks->front() == info::mpi::rank) {
			for (auto r = ranks->begin(), prev = r++; r != ranks->end(); prev = r++) {
				if (*prev != *r) {
					sOffset[n2i(*r)].push_back(goffset);
				}
			}
			goffset += _mesh->nodes->pintervals[i].end - _mesh->nodes->pintervals[i].begin;
		} else {
			rOffset[n2i(ranks->front())].push_back(0);
		}
	}

	if (!Communication::receiveLowerKnownSize(sOffset, rOffset, _mesh->neighbours)) {
		eslog::error("ESPRESO internal error: receive global offset of node intervals.\n");
	}

	_mesh->nodes->uniqueSize = goffset;
	std::vector<esint> goffsets(_mesh->neighbours.size());
	std::vector<esint> uniqueNodeOffsets = _mesh->nodes->gatherUniqueNodeDistribution();
	_mesh->nodes->uniqueOffset = uniqueNodeOffsets[info::mpi::rank];
	_mesh->nodes->uniqueTotalSize = uniqueNodeOffsets.back();

	goffset = _mesh->nodes->uniqueOffset;
	std::vector<esint> neighDistribution({ 0 });
	ranks = _mesh->nodes->iranks->cbegin();
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i, ++ranks) {
		if (_mesh->nodes->pintervals[i].sourceProcess == info::mpi::rank) {
			_mesh->nodes->pintervals[i].globalOffset = goffset;
			goffset += _mesh->nodes->pintervals[i].end - _mesh->nodes->pintervals[i].begin;
		} else {
			int nindex = n2i(_mesh->nodes->pintervals[i].sourceProcess);
			_mesh->nodes->pintervals[i].globalOffset = rOffset[nindex][goffsets[nindex]++];
			_mesh->nodes->pintervals[i].globalOffset += uniqueNodeOffsets[_mesh->nodes->pintervals[i].sourceProcess];
		}
	}
	_mesh->nodes->permute(finalpermutation);

	_mesh->nodes->dcenter.resize(_mesh->elements->ndomains);
	std::vector<Point> centers(_mesh->nodes->pintervals.size()), mins(_mesh->nodes->pintervals.size()), maxs(_mesh->nodes->pintervals.size());

	#pragma omp parallel for
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i) {
		const auto &coordinates = _mesh->nodes->coordinates->datatarray();
		Point center, min = coordinates[_mesh->nodes->pintervals[i].begin], max = coordinates[_mesh->nodes->pintervals[i].begin];
		for (esint n = _mesh->nodes->pintervals[i].begin; n < _mesh->nodes->pintervals[i].end; ++n) {
			center += coordinates[n];
			min.x = std::min(min.x, coordinates[n].x);
			min.y = std::min(min.y, coordinates[n].y);
			min.z = std::min(min.z, coordinates[n].z);
			max.x = std::max(max.x, coordinates[n].x);
			max.y = std::max(max.y, coordinates[n].y);
			max.z = std::max(max.z, coordinates[n].z);
		}
		centers[i] = center;
		mins[i] = min;
		maxs[i] = max;
	}

	_mesh->nodes->lmin = mins.front();
	_mesh->nodes->lmax = maxs.front();
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i) {
		_mesh->nodes->lmin.x = std::min(_mesh->nodes->lmin.x, mins[i].x);
		_mesh->nodes->lmin.y = std::min(_mesh->nodes->lmin.y, mins[i].y);
		_mesh->nodes->lmin.z = std::min(_mesh->nodes->lmin.z, mins[i].z);
		_mesh->nodes->lmax.x = std::max(_mesh->nodes->lmax.x, maxs[i].x);
		_mesh->nodes->lmax.y = std::max(_mesh->nodes->lmax.y, maxs[i].y);
		_mesh->nodes->lmax.z = std::max(_mesh->nodes->lmax.z, maxs[i].z);
	}

	#pragma omp parallel for
	for (size_t d = 0; d < _mesh->nodes->dintervals.size(); ++d) {
		Point center;
		esint size = 0;
		for (size_t i = 0; i < _mesh->nodes->dintervals[d].size(); ++i) {
			center += centers[_mesh->nodes->dintervals[d][i].pindex];
			size += _mesh->nodes->dintervals[d][i].end - _mesh->nodes->dintervals[d][i].begin;
		}
		_mesh->nodes->dcenter[d] = center / size;
	}

	for (size_t i = 0; i < centers.size(); i++) {
		_mesh->nodes->center += centers[i];
	}
	_mesh->nodes->center /= _mesh->nodes->size;

	MPI_Allreduce(&_mesh->nodes->lmin.x, &_mesh->nodes->min.x, 1, MPI_DOUBLE, MPI_MIN, info::mpi::comm);
	MPI_Allreduce(&_mesh->nodes->lmin.y, &_mesh->nodes->min.y, 1, MPI_DOUBLE, MPI_MIN, info::mpi::comm);
	MPI_Allreduce(&_mesh->nodes->lmin.z, &_mesh->nodes->min.z, 1, MPI_DOUBLE, MPI_MIN, info::mpi::comm);
	MPI_Allreduce(&_mesh->nodes->lmax.x, &_mesh->nodes->max.x, 1, MPI_DOUBLE, MPI_MAX, info::mpi::comm);
	MPI_Allreduce(&_mesh->nodes->lmax.y, &_mesh->nodes->max.y, 1, MPI_DOUBLE, MPI_MAX, info::mpi::comm);
	MPI_Allreduce(&_mesh->nodes->lmax.z, &_mesh->nodes->max.z, 1, MPI_DOUBLE, MPI_MAX, info::mpi::comm);

	std::vector<esint> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (esint i, esint j) { return finalpermutation[i] < finalpermutation[j]; });

	auto localremap = [&] (serializededata<esint, esint>* data) {
		if (data == NULL) {
			return;
		}
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (auto e = data->begin(t); e != data->end(t); ++e) {
				for (auto n = e->begin(); n != e->end(); ++n) {
					*n = backpermutation[*n];
				}
			}
		}
	};

	localremap(_mesh->elements->procNodes);

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->procNodes != NULL) {
			localremap(_mesh->boundaryRegions[r]->procNodes);
		}
		if (!StringCompare::caseInsensitiveEq(_mesh->boundaryRegions[r]->name, "ALL_NODES")) {
			if (_mesh->boundaryRegions[r]->nodes != NULL) {
				localremap(_mesh->boundaryRegions[r]->nodes);
				std::sort(_mesh->boundaryRegions[r]->nodes->datatarray().begin(), _mesh->boundaryRegions[r]->nodes->datatarray().end());
			}
		}
	}

	eslog::endln("MESH: NODES ARRANGED");
	eslog::checkpointln("MESH: NODES ARRANGED");
}

void MeshPreprocessing::arrangeElements()
{
	eslog::startln("MESH: ARRANGE ELEMENTS", "ARRANGE ELEMENTS");

	std::vector<esint> permutation(_mesh->elements->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	arrangeElementsPermutation(permutation);
	permuteElements(permutation, _mesh->elements->distribution);

	eslog::endln("MESH: ELEMENTS ARRANGED");
	eslog::checkpointln("MESH: ELEMENTS ARRANGED");
}

void MeshPreprocessing::arrangeElementsPermutation(std::vector<esint> &permutation)
{
	eslog::startln("MESH: ARRANGE ELEMENTS PERMUTATION", "ARRANGE ELEMENTS PERMUTATION");

	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (esint d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; ++d) {
			std::sort(
					permutation.begin() + _mesh->elements->elementsDistribution[d],
					permutation.begin() + _mesh->elements->elementsDistribution[d + 1],
					[&] (esint i, esint j) {
				return _mesh->elements->epointers->datatarray()[i]->code < _mesh->elements->epointers->datatarray()[j]->code;
			});
		}
	}

	std::vector<std::vector<esint> > iboundaries(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (esint d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; ++d) {
			if (d) {
				iboundaries[t].push_back(_mesh->elements->elementsDistribution[d]);
			}
			for (esint e = _mesh->elements->elementsDistribution[d] + 1; e < _mesh->elements->elementsDistribution[d + 1]; ++e) {
				if (_mesh->elements->epointers->datatarray()[permutation[e]]->code != _mesh->elements->epointers->datatarray()[permutation[e - 1]]->code) {
					iboundaries[t].push_back(e);
				}
			}
		}
	}
	utils::mergeThreadedUniqueData(iboundaries);

	_mesh->elements->eintervals.push_back(ElementsInterval(0, 0));
	_mesh->elements->eintervals.back().domain = _mesh->elements->firstDomain;
	_mesh->elements->eintervals.back().code = static_cast<int>(_mesh->elements->epointers->datatarray()[permutation[0]]->code);
	_mesh->elements->eintervalsDistribution.push_back(0);
	for (size_t i = 0; i < iboundaries[0].size(); i++) {
		_mesh->elements->eintervals.back().end = iboundaries[0][i];
		_mesh->elements->eintervals.push_back(ElementsInterval(iboundaries[0][i], iboundaries[0][i]));
		const std::vector<esint> &edist = _mesh->elements->elementsDistribution;
		_mesh->elements->eintervals.back().domain = std::lower_bound(edist.begin(), edist.end(), _mesh->elements->eintervals.back().begin + 1) - edist.begin() - 1 + _mesh->elements->firstDomain;
		_mesh->elements->eintervals.back().code = static_cast<int>(_mesh->elements->epointers->datatarray()[permutation[_mesh->elements->eintervals.back().begin]]->code);
		if ((_mesh->elements->eintervals.end() - 1)->domain != (_mesh->elements->eintervals.end() - 2)->domain) {
			_mesh->elements->eintervalsDistribution.push_back(_mesh->elements->eintervals.size() - 1);
		}
	}
	_mesh->elements->eintervals.back().end = _mesh->elements->size;
	_mesh->elements->eintervalsDistribution.push_back(_mesh->elements->eintervals.size());

	int elementstypes = static_cast<int>(Element::CODE::SIZE);
	if (elementstypes > 32) {
		eslog::error("ESPRESO internal error: increase elements-types synchronization buffer.\n");
	}

	int codes = 0;
	for (size_t i = 0; i < _mesh->elements->eintervals.size(); ++i) {
		codes |= 1 << _mesh->elements->eintervals[i].code;
		_mesh->elements->ecounters[_mesh->elements->eintervals[i].code] += _mesh->elements->eintervals[i].end - _mesh->elements->eintervals[i].begin;
	}

	int allcodes = 0;
	MPI_Allreduce(&codes, &allcodes, 1, MPI_INT, MPI_BOR, info::mpi::comm);

	for (int i = 0, bitmask = 1; i < elementstypes; i++, bitmask = bitmask << 1) {
		if (allcodes & bitmask) {
			_mesh->elements->ecounters[i] = Communication::exscan(_mesh->elements->ecounters[i]);
		}
	}

	eslog::endln("MESH: ELEMENTS PERMUTATION ARRANGED");
	eslog::checkpointln("MESH: ELEMENTS PERMUTATION ARRANGED");
}


void MeshPreprocessing::arrangeRegions()
{
	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	if (_mesh->elements->regions == NULL) {
		fillRegionMask();
	}

	eslog::startln("MESH: ARRANGE REGIONS", "ARRANGE REGIONS");

	size_t threads = info::env::OMP_NUM_THREADS;

	for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
		const auto &elements = _mesh->elementsRegions[r]->elements->datatarray();

		_mesh->elementsRegions[r]->eintervals = _mesh->elements->eintervals;
		for (size_t i = 0; i < _mesh->elementsRegions[r]->eintervals.size(); ++i) {
			_mesh->elementsRegions[r]->eintervals[i].begin = std::lower_bound(elements.begin(), elements.end(), _mesh->elementsRegions[r]->eintervals[i].begin) - elements.begin();
			_mesh->elementsRegions[r]->eintervals[i].end = std::lower_bound(elements.begin(), elements.end(), _mesh->elementsRegions[r]->eintervals[i].end) - elements.begin();
		}
	}

	for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
		std::vector<std::vector<esint> > unique(threads);
		_mesh->elementsRegions[r]->ueintervals = _mesh->elementsRegions[r]->eintervals;

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			const auto &regions = _mesh->elements->regions->datatarray();
			int maskSize = _mesh->elements->regionMaskSize;
			esint maskOffset = r / (8 * sizeof(esint));
			esint bit = 1 << (r % (8 * sizeof(esint)));
			std::vector<esint> mask(maskSize);
			mask[maskOffset] = bit;

			for (esint d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
				for (esint i = _mesh->elements->eintervalsDistribution[d]; i < _mesh->elements->eintervalsDistribution[d + 1]; i++) {
					size_t usize = unique[t].size();
					for (esint e = _mesh->elementsRegions[r]->eintervals[i].begin; e < _mesh->elementsRegions[r]->eintervals[i].end; ++e) {
						if (memcmp(regions.data() + e * maskSize, mask.data(), sizeof(esint) * maskSize) == 0) {
							unique[t].push_back(e);
						}
					}
					_mesh->elementsRegions[r]->ueintervals[i].begin = 0;
					_mesh->elementsRegions[r]->ueintervals[i].end = unique[t].size() - usize;
				}
			}
		}

		for (size_t i = 1; i < _mesh->elementsRegions[r]->ueintervals.size(); ++i) {
			_mesh->elementsRegions[r]->ueintervals[i].begin += _mesh->elementsRegions[r]->ueintervals[i - 1].end;
			_mesh->elementsRegions[r]->ueintervals[i].end   += _mesh->elementsRegions[r]->ueintervals[i - 1].end;
		}

		if (_mesh->elementsRegions[r]->eintervals == _mesh->elementsRegions[r]->ueintervals) {
			_mesh->elementsRegions[r]->uniqueElements = _mesh->elementsRegions[r]->elements;
		} else {
			_mesh->elementsRegions[r]->uniqueElements = new serializededata<esint, esint>(1, unique);
		}

		auto computeCountersAndNodes = [&] (const std::vector<ElementsInterval> &eintervals, const tarray<esint> &elements) {
			for (size_t i = 0; i < eintervals.size(); ++i) {
				_mesh->elementsRegions[r]->ecounters[eintervals[i].code] += eintervals[i].end - eintervals[i].begin;
			}

			std::vector<std::vector<esint> > nodes(threads);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (esint d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
					for (esint i = _mesh->elements->eintervalsDistribution[d]; i < _mesh->elements->eintervalsDistribution[d + 1]; i++) {
						if (eintervals[i].end - eintervals[i].begin > 0) {
							if (eintervals[i].end - eintervals[i].begin == _mesh->elements->eintervals[i].end - _mesh->elements->eintervals[i].begin) {
								nodes[t].insert(
										nodes[t].end(),
										(_mesh->elements->procNodes->cbegin() + _mesh->elements->eintervals[i].begin)->begin(),
										(_mesh->elements->procNodes->cbegin() + _mesh->elements->eintervals[i].end)->begin());
							} else {
								auto enodes = _mesh->elements->procNodes->cbegin() + _mesh->elements->eintervals[i].begin;
								esint prev = _mesh->elements->eintervals[i].begin;
								for (esint e = eintervals[i].begin; e < eintervals[i].end; prev = elements[e++]) {
									enodes += elements[e] - prev;
									nodes[t].insert(nodes[t].end(), enodes->begin(), enodes->end());
								}
							}
						}
					}
				}
				utils::sortAndRemoveDuplicity(nodes[t]);
			}
			utils::mergeThreadedUniqueData(nodes);
			nodes.resize(1);
			nodes.resize(threads);
			serializededata<esint, esint>::balance(1, nodes);

			_mesh->elementsRegions[r]->nodes = new serializededata<esint, esint>(1, nodes);
		};

		if (StringCompare::caseInsensitiveEq(_mesh->elementsRegions[r]->name, "ALL_ELEMENTS")) {
			computeCountersAndNodes(_mesh->elementsRegions[r]->ueintervals, _mesh->elementsRegions[r]->uniqueElements->datatarray());
		} else {
			computeCountersAndNodes(_mesh->elementsRegions[r]->eintervals, _mesh->elementsRegions[r]->elements->datatarray());
		}

		ElementsRegionStore *store = _mesh->elementsRegions[r];

		synchronizeRegionNodes(store->name, store->nodes, store->nintervals);
		computeIntervalOffsets(store->nintervals, store->uniqueOffset, store->uniqueSize, store->uniqueTotalSize);
	}

	{ // compute eregion intersection
//		std::vector<esint> nonunique(_mesh->eregion("ALL_ELEMENTS")->elements->structures() - _mesh->eregion("ALL_ELEMENTS")->uniqueElements->structures());
//		std::set_difference(
//				_mesh->eregion("ALL_ELEMENTS")->elements->datatarray().cbegin(), _mesh->eregion("ALL_ELEMENTS")->elements->datatarray().cend(),
//				_mesh->eregion("ALL_ELEMENTS")->uniqueElements->datatarray().cbegin(), _mesh->eregion("ALL_ELEMENTS")->uniqueElements->datatarray().cend(),
//				nonunique.begin());
//
//		int maskSize = _mesh->elements->regionMaskSize;
//		std::vector<std::vector<int> > tintersections(threads);
//		const auto &regions = _mesh->elements->regions->datatarray();
//		std::vector<size_t> nudistribution = tarray<esint>::distribute(threads, nonunique.size());
//
//		auto wasFound = [&] (std::vector<int> &found, const int *mask) {
//			for (size_t i = 0; i < found.size(); i += maskSize) {
//				if (memcmp(found.data() + i, mask, sizeof(int) * maskSize) == 0) {
//					return true;
//				}
//			}
//			return false;
//		};
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			for (esint e = nudistribution[t]; e < nudistribution[t + 1]; ++e) {
//				if (!wasFound(tintersections[t], regions.data() + e * maskSize)) {
//					tintersections[t].insert(tintersections[t].end(), regions.begin() + e * maskSize, regions.begin() + (e + 1) * maskSize);
//				}
//			}
//		}
//
//		std::vector<int> cintersections;
//		for (size_t t = 0; t < threads; t++) {
//			for (size_t i = 0; i < tintersections[t].size(); i += maskSize) {
//				if (!wasFound(cintersections, tintersections[t].data() + i)) {
//					cintersections.insert(cintersections.end(), tintersections[t].begin() + i, tintersections[t].begin() + i + maskSize);
//				}
//			}
//		}
//
//		std::vector<int> gintersections;
//		Communication::gatherUnknownSize(cintersections, gintersections);
//		gintersections.swap(cintersections);
//		for (size_t i = 0; i < cintersections.size(); i += maskSize) {
//			if (!wasFound(gintersections, cintersections.data() + i)) {
//				gintersections.insert(gintersections.end(), cintersections.begin() + i, cintersections.begin() + i + maskSize);
//			}
//		}
//		Communication::broadcastUnknownSize(gintersections);
//
//		std::vector<std::vector<esint> > unique(threads);
//		for (size_t i = 0; i < gintersections.size(); i += maskSize) {
//			std::string name;
//			for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
//				int maskOffset = r / (8 * sizeof(int));
//				int bit = 1 << (r % (8 * sizeof(int)));
//				if (gintersections[i + maskOffset] & bit) {
//					name += "_" + _mesh->elementsRegions[r]->name;
//				}
//			}
//			_mesh->elementsRegions.push_back(new ElementsRegionStore(name));
//
//			_mesh->elementsRegions.back()->ueintervals = _mesh->elements->eintervals;
//
//			#pragma omp parallel for
//			for (size_t t = 0; t < threads; t++) {
//				unique[t].clear();
//				for (esint d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
//					esint begin = std::lower_bound(nonunique.begin(), nonunique.end(), _mesh->elements->elementsDistribution[d]) - nonunique.begin();
//					esint end = std::lower_bound(nonunique.begin(), nonunique.end(), _mesh->elements->elementsDistribution[d + 1]) - nonunique.begin();
//					size_t usize = unique[t].size();
//					for (esint e = begin; e < end; ++e) {
//						if (memcmp(regions.data() + e * maskSize, gintersections.data() + i, sizeof(int) * maskSize) == 0) {
//							unique[t].push_back(e);
//						}
//					}
//					_mesh->elementsRegions.back()->ueintervals[d].begin = 0;
//					_mesh->elementsRegions.back()->ueintervals[d].end = unique[t].size() - usize;
//				}
//			}
//
//			_mesh->elementsRegions.back()->uniqueElements = new serializededata<esint, esint>(1, unique);
//			for (size_t i = 1; i < _mesh->elementsRegions.back()->ueintervals.size(); ++i) {
//				_mesh->elementsRegions.back()->ueintervals[i].begin += _mesh->elementsRegions.back()->ueintervals[i - 1].end;
//				_mesh->elementsRegions.back()->ueintervals[i].end   += _mesh->elementsRegions.back()->ueintervals[i - 1].end;
//			}
//		}
//
//		#pragma omp parallel for
//		for (size_t t = 0; t < threads; t++) {
//			for (esint e = nudistribution[t]; e < nudistribution[t + 1]; ++e) {
//				if (!wasFound(tintersections[t], regions.data() + e * maskSize)) {
//					tintersections[t].insert(tintersections[t].end(), regions.begin() + e * maskSize, regions.begin() + (e + 1) * maskSize);
//				}
//			}
//		}
	}

	for (size_t i = 0; i < _mesh->elements->ecounters.size(); i++) {
		if (_mesh->elements->ecounters[i]) {
			for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
				_mesh->elementsRegions[r]->ecounters[i] = Communication::exscan(_mesh->elementsRegions[r]->ecounters[i]);
			}
		}
	}

	if (_mesh->eregion("ALL_ELEMENTS")->uniqueTotalSize) {
		ElementsRegionStore* nameless = _mesh->eregion("ALL_ELEMENTS");
		_mesh->elementsRegions.push_back(new ElementsRegionStore("NAMELESS_ELEMENT_SET"));
		_mesh->elementsRegions.back()->ecounters = nameless->ecounters;
		_mesh->elementsRegions.back()->eintervals = nameless->ueintervals;
		_mesh->elementsRegions.back()->elements = new serializededata<esint, esint>(*nameless->uniqueElements);
		_mesh->elementsRegions.back()->uniqueElements = _mesh->elementsRegions.back()->elements;
		_mesh->elementsRegions.back()->nodes = new serializededata<esint, esint>(*nameless->nodes);
		_mesh->elementsRegions.back()->uniqueOffset = nameless->uniqueOffset;
		_mesh->elementsRegions.back()->uniqueSize = nameless->uniqueSize;
		_mesh->elementsRegions.back()->uniqueTotalSize = nameless->uniqueTotalSize;
		_mesh->elementsRegions.back()->nintervals = nameless->nintervals;
	}

	esint eoffset = _mesh->elements->gatherElementsProcDistribution()[info::mpi::rank];
	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->nodes == NULL) {
			std::vector<std::vector<esint> > nodes(threads);
			nodes[0] = std::vector<esint>(_mesh->boundaryRegions[r]->procNodes->datatarray().begin(), _mesh->boundaryRegions[r]->procNodes->datatarray().end());
			utils::sortAndRemoveDuplicity(nodes[0]);
			serializededata<esint, esint>::balance(1, nodes);
			_mesh->boundaryRegions[r]->nodes = new serializededata<esint, esint>(1, nodes);
		}
		BoundaryRegionStore *store = _mesh->boundaryRegions[r];
		if (!StringCompare::caseInsensitiveEq(store->name, "ALL_NODES")) {
			synchronizeRegionNodes(store->name, store->nodes, store->nintervals);
			computeIntervalOffsets(store->nintervals, store->uniqueOffset, store->uniqueSize, store->uniqueTotalSize);
		} else {
			store->nintervals = _mesh->nodes->pintervals;
			store->uniqueOffset = _mesh->nodes->uniqueOffset;
			store->uniqueSize = _mesh->nodes->uniqueSize;
			store->uniqueTotalSize = _mesh->nodes->uniqueTotalSize;
		}

		if (store->dimension == 0) {
			store->procNodes = new serializededata<esint, esint>(1, tarray<esint>(threads, 0));
			store->epointers = new serializededata<esint, Element*>(1, tarray<Element*>(threads, 0));
		} else {
			std::vector<size_t> distribution = tarray<size_t>::distribute(threads, store->procNodes->structures());
			std::vector<esint> &eDomainDistribution = _mesh->elements->elementsDistribution;
			std::vector<esint> emembership(distribution.back()), edomain(distribution.back());

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto enodes = store->procNodes->cbegin() + distribution[t];
				std::vector<esint> nlinks;
				size_t counter;
				for (size_t e = distribution[t]; e < distribution[t + 1]; ++e, ++enodes) {
					nlinks.clear();
					for (auto n = enodes->begin(); n != enodes->end(); ++n) {
						auto links = _mesh->nodes->elements->cbegin() + *n;
						nlinks.insert(nlinks.end(), links->begin(), links->end());
					}
					std::sort(nlinks.begin(), nlinks.end());
					counter = 1;
					for (size_t i = 1; i < nlinks.size(); ++i) {
						if (nlinks[i - 1] == nlinks[i]) {
							++counter;
							if (counter == enodes->size() && eoffset <= nlinks[i]) {
								emembership[e] = nlinks[i];
								edomain[e] = std::lower_bound(eDomainDistribution.begin(), eDomainDistribution.end(), nlinks[i] + 1 - eoffset) - eDomainDistribution.begin() - 1;
								break;
							}
						} else {
							counter = 1;
						}
					}
				}
			}

			std::vector<esint> permutation(store->procNodes->structures());
			std::iota(permutation.begin(), permutation.end(), 0);
			std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
				if (edomain[i] == edomain[j]) {
					if (store->epointers->datatarray()[i] == store->epointers->datatarray()[j]) {
						return emembership[i] < emembership[j];
					}
					return store->epointers->datatarray()[i] < store->epointers->datatarray()[j];
				}
				return edomain[i] < edomain[j];
			});

			std::vector<size_t> edistribution;
			for (auto i = _mesh->elements->elementsDistribution.begin(); i != _mesh->elements->elementsDistribution.end(); ++i) {
				auto it = std::lower_bound(permutation.begin(), permutation.end(), *i, [&] (esint i, esint d) { return emembership[i] - eoffset < d; });
				edistribution.push_back(it - permutation.begin());
			}

			std::vector<size_t> tdistribution;
			for (size_t t = 0; t < _mesh->elements->domainDistribution.size(); t++) {
				tdistribution.push_back(edistribution[_mesh->elements->domainDistribution[t]]);
			}

			store->permute(permutation, tdistribution);

			std::vector<std::vector<esint> > iboundaries(threads);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (esint d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
					iboundaries[t].push_back(edistribution[d]);
					for (size_t e = edistribution[d] + 1; e < edistribution[d + 1]; ++e) {
						if (store->epointers->datatarray()[e]->code != store->epointers->datatarray()[e - 1]->code) {
							iboundaries[t].push_back(e);
						}
					}
				}
			}
			iboundaries.back().push_back(edistribution.back());
			utils::mergeThreadedUniqueData(iboundaries);

			store->eintervalsDistribution.push_back(0);
			esint lastDomain = 0;
			for (size_t i = 0; i < iboundaries[0].size() - 1; i++) {
				store->eintervals.push_back(ElementsInterval(iboundaries[0][i], iboundaries[0][i + 1]));
				store->eintervals.back().code = static_cast<int>(store->epointers->datatarray()[iboundaries[0][i]]->code);
				store->eintervals.back().domain = edomain[permutation[iboundaries[0][i]]];
				if (store->eintervals.back().domain != lastDomain) {
					store->eintervalsDistribution.insert(
							store->eintervalsDistribution.end(),
							store->eintervals.back().domain - lastDomain,
							store->eintervals.size() - 1);
				}
				lastDomain = store->eintervals.back().domain;
			}
			store->eintervalsDistribution.insert(
					store->eintervalsDistribution.end(),
					_mesh->elements->ndomains - lastDomain,
					store->eintervals.size());

			int codes = 0;
			for (size_t i = 0; i < store->eintervals.size(); ++i) {
				store->ecounters[store->eintervals[i].code] += store->eintervals[i].end - store->eintervals[i].begin;
				codes |= 1 << store->eintervals[i].code;
			}

			int allcodes = 0;
			MPI_Allreduce(&codes, &allcodes, 1, MPI_INT, MPI_BOR, info::mpi::comm);

			for (size_t i = 0, bitmask = 1; i < _mesh->elements->ecounters.size(); i++, bitmask = bitmask << 1) {
				if (allcodes & bitmask) {
					store->ecounters[i] = Communication::exscan(store->ecounters[i]);
				}
			}
		}
	}

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->dimension) {
			computeRegionArea(_mesh->boundaryRegions[r]);
		}
	}

	eslog::endln("MESH: REGIONS ARRANGED");
}

void MeshPreprocessing::fillRegionMask()
{
	eslog::startln("MESH: FILL REGION MASK", "FILL REGIONS MASKS");

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > eregions(threads);

	// regions are transfered via mask
	int regionsBitMaskSize = _mesh->elementsRegions.size() / (8 * sizeof(esint)) + (_mesh->elementsRegions.size() % (8 * sizeof(esint)) ? 1 : 0);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint maskOffset = 0;
		eregions[t].resize(regionsBitMaskSize * (_mesh->elements->distribution[t + 1] - _mesh->elements->distribution[t]));
		for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(esint));
			esint bit = 1 << (r % (8 * sizeof(esint)));

			const auto &elements = _mesh->elementsRegions[r]->elements->datatarray();
			auto begin = std::lower_bound(elements.begin(), elements.end(), _mesh->elements->distribution[t]);
			auto end = std::lower_bound(elements.begin(), elements.end(), _mesh->elements->distribution[t + 1]);
			for (auto i = begin; i != end; ++i) {
				eregions[t][(*i - _mesh->elements->distribution[t]) * regionsBitMaskSize + maskOffset] |= bit;
			}
		}
	}

	_mesh->elements->regionMaskSize = regionsBitMaskSize;
	_mesh->elements->regions = new serializededata<esint, esint>(regionsBitMaskSize, eregions);

	eslog::endln("MESH: REGION MASK FILLED");
}

void MeshPreprocessing::synchronizeRegionNodes(const std::string &name, serializededata<esint, esint>* &rnodes, std::vector<ProcessInterval> &nintervals)
{
	eslog::start("MESH: SYNCHRONIZE REGION", "REGION SYNCHRONIZATION");
	eslog::param("REGION", name.c_str());
	eslog::ln();

	const auto &nodes = rnodes->datatarray();

	auto n2i = [ & ] (int neighbour) -> int {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbour) - _mesh->neighbours.begin();
	};

	std::vector<std::vector<esint> > sBuffer(_mesh->neighbours.size()), rBuffer(_mesh->neighbours.size());

	int prevRank;
	auto iranks = _mesh->nodes->iranks->cbegin();
	nintervals = _mesh->nodes->pintervals;
	for (size_t i = 0; i < nintervals.size(); ++i, ++iranks) {
		nintervals[i].begin = std::lower_bound(nodes.begin(), nodes.end(), nintervals[i].begin) - nodes.begin();
		nintervals[i].end = std::lower_bound(nodes.begin(), nodes.end(), nintervals[i].end) - nodes.begin();

		prevRank = -1;
		esint isize = nintervals[i].end - nintervals[i].begin;
		for (auto rank = iranks->begin(); rank != iranks->end(); ++rank) {
			if (*rank != info::mpi::rank && *rank != prevRank) {
				sBuffer[n2i(*rank)].push_back(_mesh->nodes->pintervals[i].globalOffset);
				sBuffer[n2i(*rank)].push_back(isize);
				for (auto n = nodes.begin() + nintervals[i].begin; n != nodes.begin() + nintervals[i].end; ++n) {
					sBuffer[n2i(*rank)].push_back(*n - _mesh->nodes->pintervals[i].begin);
				}
				prevRank = *rank;
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, _mesh->neighbours)) {
		eslog::error("ESPRESO internal error: exchange element region nodes.\n");
	}

	std::vector<esint> nnodes;
	for (size_t neigh = 0; neigh < _mesh->neighbours.size(); ++neigh) {
		if (sBuffer[neigh] != rBuffer[neigh]) {
			auto it = rBuffer[neigh].begin();
			while (it != rBuffer[neigh].end()) {
				esint globalOffset = *it; ++it;
				esint isize = *it; ++it;
				if (isize) {
					const ProcessInterval &interval = *std::lower_bound(_mesh->nodes->pintervals.begin(), _mesh->nodes->pintervals.end(), globalOffset, [&] (const ProcessInterval &interval, esint offset) {
						return interval.globalOffset < offset;
					});
					for (esint n = 0; n < isize; ++n, ++it) {
						nnodes.push_back(interval.begin + *it);
					}
				}
			}
		}
	}
	if (nnodes.size()) {
		nnodes.insert(nnodes.end(), nodes.begin(), nodes.end());
		utils::sortAndRemoveDuplicity(nnodes);
		delete rnodes;
		rnodes = new serializededata<esint, esint>(1, nnodes);

		nintervals = _mesh->nodes->pintervals;
		for (size_t i = 0; i < nintervals.size(); ++i) {
			nintervals[i].begin = std::lower_bound(nnodes.begin(), nnodes.end(), nintervals[i].begin) - nnodes.begin();
			nintervals[i].end = std::lower_bound(nnodes.begin(), nnodes.end(), nintervals[i].end) - nnodes.begin();
		}
	}

	eslog::end("MESH: REGION SYNCHRONIZED");
	eslog::param("REGION", name.c_str());
	eslog::ln();
}

void MeshPreprocessing::computeBoundaryElementsFromNodes(BoundaryRegionStore *bregion, int elementDimension)
{
	if (_mesh->nodes->elements == NULL) {
		linkNodesAndElements();
	}

	if (_mesh->elements->neighbors == NULL) {
		computeElementsNeighbors();
	}

	eslog::start("MESH: BOUNDARY FROM NODES", "BOUNDARY FROM NODES");
	eslog::param("BOUNDARY", bregion->name.c_str());
	eslog::ln();

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<std::pair<esint, esint> > > elements(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = bregion->nodes->begin(t)->begin(); n != bregion->nodes->end(t)->begin(); ++n) {
			auto links = _mesh->nodes->elements->begin() + *n;
			for (auto e = links->begin(); e != links->end(); ++e) {
				elements[t].push_back(std::make_pair(*e, *n));
			}
		}
	}

	std::vector<size_t> distribution = { 0, elements[0].size() };
	for (size_t t = 1; t < threads; t++) {
		elements[0].insert(elements[0].end(), elements[t].begin(), elements[t].end());
		distribution.push_back(elements[0].size());
	}

	utils::sortWithInplaceMerge(elements[0], distribution);

	std::vector<esint> edistribution = _mesh->elements->gatherElementsProcDistribution();
	esint ebegin = edistribution[info::mpi::rank];
	esint eend = edistribution[info::mpi::rank + 1];

	auto begin = std::lower_bound(elements[0].begin(), elements[0].end(), ebegin,
			[] (const std::pair<esint, esint> &p, esint e) { return p.first < e; });
	auto end = std::lower_bound(elements[0].begin(), elements[0].end(), eend,
			[] (const std::pair<esint, esint> &p, esint e) { return p.first < e; });

	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, end - begin);
	for (size_t t = 1; t < threads; t++) {
		while (
				begin + tdistribution[t] < end && begin <= begin + tdistribution[t] - 1 &&
				(begin + tdistribution[t] - 1)->first == (begin + tdistribution[t])->first) {

			++tdistribution[t];
		}
	}

	std::vector<std::vector<esint> > edist(threads), edata(threads), ecode(threads);
	std::vector<std::vector<Element*> > epointers(threads);

	int rsize = _mesh->elements->regionMaskSize;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> nodes, facenodes, lowerElements, lenodes;
		std::vector<esint> tdist, tdata, tcode;
		if (t == 0) {
			tdist.push_back(0);
		}

		int nface;
		esint element, neighbor, prev = 0;
		auto enodes = _mesh->elements->procNodes->cbegin();
		auto neighbors = _mesh->elements->neighbors->cbegin();
		const auto &regions = _mesh->elements->regions->datatarray();

		for (size_t e = tdistribution[t]; e < tdistribution[t + 1]; e++) {
			nodes.push_back((begin + e)->second);
			if ((e + 1 == tdistribution[t + 1] || (begin + e + 1)->first != (begin + e)->first)) {

				element = (begin + e)->first - ebegin;
				utils::sortAndRemoveDuplicity(nodes);

				enodes += element - prev;
				neighbors += element - prev;
				prev = element;

				const auto &fpointers = _mesh->elements->epointers->datatarray()[element]->facepointers->datatarray();
				auto fnodes = _mesh->elements->epointers->datatarray()[element]->faces->cbegin();
				nface = 0;
				for (auto f = fpointers.begin(); f != fpointers.end(); ++f, ++fnodes, ++nface) {

					auto addFace = [&] () {
						for (auto n = fnodes->begin(); n != fnodes->end(); ++n) {
							tdata.push_back(enodes->at(*n));
						}
						tdist.push_back(tdata.size());
						tcode.push_back((esint)(*f)->code);
					};

					if ((int)nodes.size() >= (*f)->nodes) {
						for (auto n = fnodes->begin(); n != fnodes->end(); ++n) {
							facenodes.push_back(enodes->at(*n));
						}
						std::sort(facenodes.begin(), facenodes.end());
						if (std::includes(nodes.begin(), nodes.end(), facenodes.begin(), facenodes.end())) {
							neighbor = neighbors->at(nface);
							if (neighbor == -1) {
								addFace();
							} else if (element + ebegin < neighbor) {
								if (ebegin <= neighbor && neighbor < eend) {
									neighbor -= ebegin;
									if (memcmp(regions.data() + element * rsize, regions.data() + neighbor * rsize, sizeof(esint) * rsize) != 0) {
										addFace();
									}
								} else {
									neighbor = std::lower_bound(_mesh->halo->IDs->datatarray().begin(), _mesh->halo->IDs->datatarray().end(), neighbor) - _mesh->halo->IDs->datatarray().begin();
									if (memcmp(regions.data() + element * rsize, _mesh->halo->regions->datatarray().data() + neighbor * rsize, sizeof(esint) * rsize) != 0) {
										addFace();
									}
								}
							}
						}
						facenodes.clear();
					}
				}
				nodes.clear();
			}
		}

		edist[t].swap(tdist);
		edata[t].swap(tdata);
		ecode[t].swap(tcode);
	}

	utils::threadDistributionToFullDistribution(edist);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = 0; e < ecode[t].size(); e++) {
			epointers[t].push_back(Mesh::edata + ecode[t][e]);
		}
	}

	bregion->procNodes = new serializededata<esint, esint>(edist, edata);
	bregion->epointers = new serializededata<esint, Element*>(1, epointers);
	bregion->dimension = elementDimension;

	eslog::end("MESH: BOUNDARY FROM NODES COMPUTED");
	eslog::param("BOUNDARY", bregion->name.c_str());
	eslog::ln();
}

void MeshPreprocessing::computeIntervalOffsets(std::vector<ProcessInterval> &intervals, esint &uniqueOffset, esint &uniqueSize, esint &uniqueTotalSize)
{
	auto n2i = [ & ] (int neighbour) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbour) - _mesh->neighbours.begin();
	};

	std::vector<std::vector<esint> > sOffset(_mesh->neighbours.size()), rOffset(_mesh->neighbours.size());
	auto ranks = _mesh->nodes->iranks->cbegin();
	esint myoffset = 0;
	for (size_t i = 0; i < intervals.size(); ++i, ++ranks) {
		intervals[i].sourceProcess = ranks->front();
		if (ranks->front() == info::mpi::rank) {
			for (auto r = ranks->begin(), prev = r++; r != ranks->end(); prev = r++) {
				if (*prev != *r) {
					sOffset[n2i(*r)].push_back(myoffset);
				}
			}
			myoffset += intervals[i].end - intervals[i].begin;
		} else {
			rOffset[n2i(ranks->front())].push_back(0);
		}
	}

	if (!Communication::receiveLowerKnownSize(sOffset, rOffset, _mesh->neighbours)) {
		eslog::error("ESPRESO internal error: receive global offset of intervals.\n");
	}

	uniqueSize = myoffset;
	std::vector<esint> uniqueOffsets = Store::gatherDistribution(uniqueSize);
	uniqueOffset = uniqueOffsets[info::mpi::rank];
	uniqueTotalSize = uniqueOffsets.back();

	std::vector<esint> roffsetsIndex(_mesh->neighbours.size());

	myoffset = uniqueOffset;
	ranks = _mesh->nodes->iranks->cbegin();
	for (size_t i = 0; i < intervals.size(); ++i, ++ranks) {
		if (intervals[i].sourceProcess == info::mpi::rank) {
			intervals[i].globalOffset = myoffset;
			myoffset += intervals[i].end - intervals[i].begin;
		} else {
			int nindex = n2i(intervals[i].sourceProcess);
			intervals[i].globalOffset = rOffset[nindex][roffsetsIndex[nindex]++];
			intervals[i].globalOffset += uniqueOffsets[intervals[i].sourceProcess];
		}
	}
}

