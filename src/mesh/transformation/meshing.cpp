
#include "transformations.h"

#include "../mesh.h"
#include "../elements/elementstore.h"
#include "../store/domainstore.h"
#include "../store/boundarystore.h"
#include "../store/regionstore.h"

#include "../../basis/point/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../../config/ecf/environment.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

void Transformation::computeElementCenters(Mesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::computation of elements centers started.";

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<Point> > centers(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		Point center;
		for (auto nodes = mesh._elems->nodes->cbegin(t); nodes != mesh._elems->nodes->cend(t); ++nodes) {
			center = Point();
			for (auto n = nodes->begin(); n != nodes->end(); ++n) {
				center += mesh._nodes->coordinates->datatarray()[*n];
			}
			center /= nodes->size();
			centers[t].push_back(center);
		}
	}

	mesh._elems->coordinates = new serializededata<eslocal, Point>(1, centers);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of elements centers finished.";
}

void Transformation::computeDomainsCenters(Mesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::computation of domain centers started.";

//	if (mesh._nodes->domains == NULL) {
//		Transformation::assignDomainsToNodes(mesh);
//	}
//
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	std::vector<std::vector<Point> > centers(threads, std::vector<Point>(mesh._domains->structures()));
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		auto domains = mesh._nodes->domains->cbegin(t);
//		for (auto n = mesh._nodes->coordinates->cbegin(t); n != mesh._nodes->coordinates->cend(t); ++n, ++domains) {
//			for (auto d = domains->begin(); d != domains->end(); ++d) {
//				;
//			}
//		}
//	}
//
//	mesh._elems->coordinates = new serializededata<eslocal, Point>(1, centers);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of domain centers finished.";
}

void Transformation::reindexNodes(Mesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::re-index nodes started.";

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<esglobal> myIDs(threads);

	auto n2i = [&] (int neighbor) {
		return std::lower_bound(mesh._neighbours.begin(), mesh._neighbours.end(), neighbor) - mesh._neighbours.begin();
	};

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esglobal myID = 0;
		for (auto ranks = mesh._nodes->ranks->cbegin(t); ranks != mesh._nodes->ranks->cend(t); ++ranks) {
			if (ranks->front() == environment->MPIrank) {
				++myID;
			}
		}
		myIDs[t] = myID;
	}

	esglobal IDoffset = Esutils::sizesToOffsets(myIDs);
	if (!Communication::exscan(IDoffset)) {
		ESINFO(ERROR) << "ESPRESO internal error: exscan IDOffsets.";
	}

	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > sIDMap(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(mesh._neighbours.size()));
	std::vector<std::vector<std::pair<esglobal, esglobal> > > IDMap(mesh._neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = mesh._nodes->ranks->cbegin(t);
		auto &IDs = mesh._nodes->IDs->datatarray();

		esglobal ID, myID = myIDs[t] + IDoffset;
		for (auto n = mesh._nodes->distribution[t]; n < mesh._nodes->distribution[t + 1]; ++n, ++ranks) {
			ID = IDs[n];
			if (ranks->front() == environment->MPIrank) {
				IDs[n] = myID++;
			}
			for (auto r = ranks->begin(); r != ranks->end(); ++r) {
				if (*r < environment->MPIrank) {
					sIDMap[t][n2i(*r)].push_back(std::make_pair(ID, (esglobal)n));
				}
				if (*r > environment->MPIrank) {
					sIDMap[t][n2i(*r)].push_back(std::make_pair(ID, myID));
				}
			}
		}
	}

	Esutils::mergeThreadedUniqueData(sIDMap);

	for (size_t n = 0; n < mesh._neighbours.size(); n++) {
		if (mesh._neighbours[n] < environment->MPIrank) {
			IDMap[n].resize(sIDMap[0][n].size());
		}
	}

	if (!Communication::receiveLowerKnownSize(sIDMap[0], IDMap, mesh._neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: receive IDs from lower ranks.";
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = mesh._nodes->ranks->cbegin(t);
		auto &IDs = mesh._nodes->IDs->datatarray();

		for (auto n = mesh._nodes->distribution[t]; n < mesh._nodes->distribution[t + 1]; ++n, ++ranks) {
			if (ranks->front() < environment->MPIrank) {
				IDs[n] = std::lower_bound(IDMap[n2i(ranks->front())].begin(), IDMap[n2i(ranks->front())].end(), std::make_pair(IDs[n], 0))->second;
			}
		}
	}

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::re-index nodes finished.";
}

void Transformation::arrangeNodes(Mesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::arrange nodes started.";

	if (mesh._nodes->domains == NULL) {
		Transformation::assignDomainsToNodes(mesh);
	}

	if (mesh._processBoundaries->nodes == NULL) {
		Transformation::computeProcessBoundaries(mesh);
	}

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<esglobal> domainBoundaries = mesh._domains->gatherElementDistribution();
	std::vector<esglobal> processBoundaries = mesh._elems->gatherElementDistrubution();

	mesh._domains->offset = std::lower_bound(domainBoundaries.begin(), domainBoundaries.end(), processBoundaries[environment->MPIrank] + 1) - domainBoundaries.begin() - 1;

	eslocal externalNodeCount = 0, boundaryNodeCount = mesh._processBoundaries->nodesIntervals.back().end;
	for (size_t i = 0; i < mesh._processBoundaries->nodesIntervals.size() && mesh._processBoundaries->nodesIntervals[i].neighbors.front() == -1; ++i) {
		externalNodeCount = mesh._processBoundaries->nodesIntervals[i].end;
	}

	std::vector<eslocal> permutation;
	permutation.reserve(mesh._nodes->size);
	permutation.insert(permutation.end(), mesh._processBoundaries->nodes->datatarray().data(), mesh._processBoundaries->nodes->datatarray().data() + boundaryNodeCount);

	std::vector<eslocal> indices(mesh._processBoundaries->nodesIntervals.size());
	for (size_t i = 0; i < mesh._processBoundaries->nodesIntervals.size(); ++i) {
		indices[i] = mesh._processBoundaries->nodesIntervals[i].begin;
	}
	for (eslocal n = 0; n < (eslocal)mesh._nodes->size; ++n) {
		for (size_t i = 0; i < indices.size(); ++i) {
			if (indices[i] < mesh._processBoundaries->nodesIntervals[i].end && permutation[indices[i]] == n) {
				indices[i]++;
				break;
			}
			if (i + 1 == indices.size()) {
				permutation.push_back(n);
			}
		}
	}

	auto comp = [&] (eslocal i, eslocal j) {
		auto di = mesh._nodes->domains->cbegin() + i;
		auto dj = mesh._nodes->domains->cbegin() + j;

		if (di->size() == dj->size()) {
			for (size_t d = 0; d < di->size(); d++) {
				if ((*di)[d] != (*dj)[d]) {
					return (*di)[d] < (*dj)[d];
				}
			}
		}
		return di->size() > dj->size();
	};

	std::sort(permutation.begin(), permutation.begin() + externalNodeCount, comp);
	std::sort(permutation.begin() + externalNodeCount, permutation.begin() + boundaryNodeCount, comp);
	std::sort(permutation.begin() + boundaryNodeCount, permutation.end(), comp);

	std::vector<EInterval> nintervals;

	auto splitinterval = [&] (eslocal boundary) {
		auto it = std::lower_bound(nintervals.begin(), nintervals.end(), boundary, [] (EInterval &internal, eslocal n) { return internal.end < n; });
		if (boundary < it->end) {
			size_t i = it - nintervals.begin();
			nintervals.insert(it, *it);
			nintervals[i].end = boundary;
			nintervals[i + 1].begin = boundary;
		}
	};

	Transformation::computeIntervals(nintervals, *mesh._nodes->domains, mesh._nodes->distribution, permutation);
	splitinterval(externalNodeCount);
	splitinterval(boundaryNodeCount);
	auto iti = nintervals.begin();
	while (iti->end <= externalNodeCount) {
		iti->neighbors.insert(iti->neighbors.begin(), -1);
		++iti;
	}

	#pragma omp parallel for
	for (size_t i = 0; i < nintervals.size(); ++i) {
		std::sort(permutation.begin() + nintervals[i].begin, permutation.begin() + nintervals[i].end, [&] (eslocal i, eslocal j) {
			return mesh._nodes->IDs->datatarray()[i] < mesh._nodes->IDs->datatarray()[j];
		});
	}

	// 0 -> external + inner, 1 -> inner, 2 -> external, 3 -> rest
	auto getBoundary = [] (const std::vector<int> &neighbors) {
		if (neighbors.front() == -1) {
			if (neighbors.size() > 2) {
				return 0;
			} else {
				return 2;
			}
		} else {
			if (neighbors.size() > 1) {
				return 1;
			} else {
				return 3;
			}
		}
	};

	std::sort(nintervals.begin(), nintervals.end(), [&] (EInterval &ei1, EInterval &ei2) {
		int b1 = getBoundary(ei1.neighbors), b2 = getBoundary(ei2.neighbors);
		if (b1 == b2) {
			if (ei1.neighbors.size() == ei2.neighbors.size()) {
				return ei1.neighbors.size() > ei2.neighbors.size();
			}
			return ei1.neighbors < ei1.neighbors;
		}
		return b1 < b2;
	});

	std::vector<eslocal> finalpermutation;
	finalpermutation.reserve(permutation.size());

	eslocal ioffset = 0;
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = mesh._domains->domainDistribution[t]; d < mesh._domains->domainDistribution[t + 1]; d++) {
			for (size_t i = 0; i < nintervals.size(); i++) {
				if (*std::lower_bound(nintervals[i].neighbors.begin(), nintervals[i].neighbors.end(), mesh._domains->offset) == mesh._domains->offset + d) {
					finalpermutation.insert(finalpermutation.end(), permutation.begin() + nintervals[i].begin, permutation.begin() + nintervals[i].end);
					eslocal isize = nintervals[i].end - nintervals[i].begin;
					nintervals[i].begin = nintervals[i].clusterOffset = ioffset;
					ioffset += isize;
					nintervals[i].end = ioffset;
				}
			}
		}
	}

	mesh.computeIntervalsOffsets(
			nintervals,
			[&] (eslocal i) { return nintervals[i].end - nintervals[i].begin; },
			[&] (eslocal i, esglobal offset) { nintervals[i].globalOffset = offset; });

	mesh._domains->nodesIntervals = nintervals;
	mesh._domains->domainNodesIntervals.resize(mesh._domains->size);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (eslocal d = mesh._domains->domainDistribution[t]; d < mesh._domains->domainDistribution[t + 1]; d++) {
			eslocal doffset = 0;
			for (size_t i = 0; i < nintervals.size(); i++) {
				if (std::binary_search(nintervals[i].neighbors.begin(), nintervals[i].neighbors.end(), mesh._domains->offset + d)) {
					mesh._domains->domainNodesIntervals[d].push_back(nintervals[i]);
					mesh._domains->domainNodesIntervals[d].back().globalDomainOffset =
							std::lower_bound(nintervals[i].neighbors.begin(), nintervals[i].neighbors.end(), mesh._domains->offset + d) - nintervals[i].neighbors.begin();
					mesh._domains->domainNodesIntervals[d].back().localDomainOffset =
							std::lower_bound(nintervals[i].neighbors.begin(), nintervals[i].neighbors.end(), mesh._domains->offset) - nintervals[i].neighbors.begin();
					mesh._domains->domainNodesIntervals[d].back().localDomainOffset = mesh._domains->domainNodesIntervals[d].back().globalDomainOffset - mesh._domains->domainNodesIntervals[d].back().localDomainOffset;
					if (nintervals[i].neighbors.front() == -1) {
						--mesh._domains->domainNodesIntervals[d].back().globalDomainOffset;
					}
					mesh._domains->domainNodesIntervals[d].back().domainOffset = doffset;
					doffset += nintervals[i].end - nintervals[i].begin;
				}
			}
		}
	}

	mesh._nodes->permute(finalpermutation);

	std::vector<eslocal> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (eslocal i, eslocal j) { return finalpermutation[i] < finalpermutation[j]; });

	auto localremap = [&] (serializededata<eslocal, eslocal>* data) {
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

	localremap(mesh._elems->nodes);
	localremap(mesh._processBoundaries->nodes);
	localremap(mesh._domainsBoundaries->nodes);

	for (size_t r = 0; r < mesh._regions.size(); r++) {
		localremap(mesh._regions[r]->nodes);
		std::sort(mesh._regions[r]->nodes->datatarray().begin(), mesh._regions[r]->nodes->datatarray().end());

		{
			auto offset = mesh._regions[r]->nodes->datatarray().cbegin();
			auto cend = mesh._regions[r]->nodes->datatarray().cend();
			for (size_t i = 0; i < mesh._domains->nodesIntervals.size(); i++) {
				auto begin = offset = std::lower_bound(offset, cend, mesh._domains->nodesIntervals[i].begin);
				auto end   = offset = std::lower_bound(offset, cend, mesh._domains->nodesIntervals[i].end);
				mesh._regions[r]->nodesIntervals.push_back(mesh._domains->nodesIntervals[i]);
				mesh._regions[r]->nodesIntervals.back().begin = begin - mesh._regions[r]->nodes->datatarray().cbegin();
				mesh._regions[r]->nodesIntervals.back().end   = end   - mesh._regions[r]->nodes->datatarray().cbegin();
			}
		}

		mesh._regions[r]->domainNodesIntervals.resize(mesh._domains->size);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (eslocal d = mesh._domains->domainDistribution[t]; d < mesh._domains->domainDistribution[t + 1]; d++) {
				auto offset = mesh._regions[r]->nodes->datatarray().cbegin();
				auto cend = mesh._regions[r]->nodes->datatarray().cend();
				for (size_t i = 0; i < mesh._domains->domainNodesIntervals[d].size(); i++) {
					auto begin = offset = std::lower_bound(offset, cend, mesh._domains->domainNodesIntervals[d][i].begin);
					auto end   = offset = std::lower_bound(offset, cend, mesh._domains->domainNodesIntervals[d][i].end);
					mesh._regions[r]->domainNodesIntervals[d].push_back(mesh._domains->domainNodesIntervals[d][i]);
					mesh._regions[r]->domainNodesIntervals[d].back().begin = begin - mesh._regions[r]->nodes->datatarray().cbegin();
					mesh._regions[r]->domainNodesIntervals[d].back().end   = end   - mesh._regions[r]->nodes->datatarray().cbegin();
				}
			}
		}

	}

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::arrange nodes finished.";
}
