
#include "transformations.h"

#include "../mesh.h"
#include "../elements/element.h"
#include "../elements/elementstore.h"
#include "../store/domainstore.h"
#include "../store/boundarystore.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../config/ecf/environment.h"

#include <numeric>
#include <algorithm>
#include <iostream>

#include "../../basis/utilities/communication.h"

using namespace espreso;

template <typename Tdual>
void Transformation::computeBoundaries(Mesh &mesh,
		serializededata<eslocal, Tdual>        *elementDual,
		esglobal                               dualOffset,
		const std::vector<esglobal>            &IDBoundaries,
		std::vector<std::vector<eslocal> >     *elementData,
		std::vector<std::vector<eslocal> >     *faceDistribution,
		std::vector<std::vector<eslocal> >     *faceData,
		std::vector<std::vector<Element*> > *faceCodes,
		std::vector<std::vector<int> >         *faceNeighbors)
{
	size_t threads = environment->OMP_NUM_THREADS;

	esglobal eoffset = mesh._elems->IDs->cbegin()->front();

	if (faceDistribution != NULL) {
		faceDistribution->front().push_back(0);
	}
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esglobal> common;
		size_t ncommons, counter;
		int me, neighbor;
		auto dual = elementDual->cbegin(t);
		auto epointer = mesh._elems->epointers->cbegin(t);
		esglobal eID = mesh._elems->distribution[t];
		auto IDpointer = std::lower_bound(IDBoundaries.begin(), IDBoundaries.end(), eID + eoffset + 1) - 1;
		esglobal begine = *IDpointer, ende = *(IDpointer + 1);
		me = IDpointer - IDBoundaries.begin();

		for (auto e = mesh._elems->nodes->cbegin(t); e != mesh._elems->nodes->cend(t); ++e, ++dual, ++epointer, ++eID) {
			if (eID + eoffset >= ende) {
				++IDpointer;
				++me;
				begine = *IDpointer;
				ende = *(IDpointer + 1);
			}
			if (dual->size() < epointer->front()->faces->structures() || dual->front() + dualOffset < begine || dual->back() + dualOffset >= ende) {

				auto facepointer = epointer->front()->facepointers->cbegin(t);
				for (auto face = epointer->front()->faces->cbegin(t); face != epointer->front()->faces->cend(t); ++face, ++facepointer) {

					neighbor = -1;
					common.clear();
					for (auto n = face->begin(); n != face->end(); ++n) {
						auto nelements = mesh._nodes->elems->cbegin() + (*e)[*n];
						for (auto ne = nelements->begin(); ne != nelements->end(); ++ne) {
							common.push_back(*ne);
						}
					}
					std::sort(common.begin(), common.end());

					ncommons = counter = 0;
					for (size_t i = 1; i < common.size(); i++) {
						if (common[i - 1] == common[i]) {
							++counter;
						} else {
							if (face->size() == counter + 1) {
								if (begine <= common[i - 1] && common[i - 1] < ende) {
									++ncommons;
								} else {
									neighbor = std::lower_bound(IDBoundaries.begin(), IDBoundaries.end(), common[i - 1] + 1) - IDBoundaries.begin() - 1;
								}
							}
							counter = 0;
						}
					}
					if (face->size() == counter + 1) {
						if (begine <= common.back() && common.back() < ende) {
							++ncommons;
						} else {
							neighbor = std::lower_bound(IDBoundaries.begin(), IDBoundaries.end(), common.back() + 1) - IDBoundaries.begin() - 1;
						}
					}

					if (ncommons == 1) {
						if (faceDistribution != NULL) {
							(*faceDistribution)[t].push_back(face->size());
							(*faceCodes)[t].push_back(facepointer->front());
							if ((*faceDistribution)[t].size() > 1) {
								(*faceDistribution)[t].back() += *((*faceDistribution)[t].end() - 2);
							}
							for (auto n = face->begin(); n != face->end(); ++n) {
								(*faceData)[t].push_back((*e)[*n]);
							}
							(*faceNeighbors)[t].push_back(neighbor);
							(*faceNeighbors)[t].push_back(me);
						}
						if (elementData != NULL) {
							if (neighbor != -1 && ((*elementData)[t].size() == 0 || (*elementData)[t].back() != eID)) {
								(*elementData)[t].push_back(eID);
							}
						}
					}
				}
			}
		}
	}
}

void Transformation::distributeElementsToIntervals(Mesh &mesh,
		BoundaryStore*                         &boundaries,
		const std::vector<esglobal>            &IDBoundaries,
		std::vector<std::vector<eslocal> >     &elementData)
{
	serializededata<eslocal, eslocal>::balance(1, elementData);
	boundaries->elems = new serializededata<eslocal, eslocal>(1, elementData);
	std::vector<eslocal> permutation(boundaries->elems->structures());
	std::iota(permutation.begin(), permutation.end(), 0);

	auto moveToNextNeigh = [&] (const eslocal* &ne, const eslocal *end) {
		while (++ne != end) {
			int neigh = std::lower_bound(IDBoundaries.begin(), IDBoundaries.end(), *ne + 1) - IDBoundaries.begin() - 1;
			if (neigh != environment->MPIrank) {
				return neigh;
			}
		}
		return -1;
	};
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		auto di = mesh._elems->dual->cbegin() + boundaries->elems->datatarray()[i];
		auto dj = mesh._elems->dual->cbegin() + boundaries->elems->datatarray()[j];

		auto diit = di->begin();
		auto djit = dj->begin();
		int ni = moveToNextNeigh(--diit, di->end());
		int nj = moveToNextNeigh(--djit, dj->end());
		while (diit != di->end() && djit != dj->end()) {
			if (ni == -1 && nj == -1) {
				break;
			}
			if (ni == -1) {
				return true;
			}
			if (nj == -1) {
				return false;
			}
			if (ni != nj) {
				return ni < nj;
			}
			ni = moveToNextNeigh(diit, di->end());
			nj = moveToNextNeigh(djit, dj->end());
		}
		return boundaries->elems->datatarray()[i] < boundaries->elems->datatarray()[j];
	});

	boundaries->elems->permute(permutation);
}

void Transformation::distributeFacesToIntervals(Mesh &mesh,
		BoundaryStore*                         &boundaries,
		std::vector<std::vector<eslocal> >     &faceDistribution,
		std::vector<std::vector<eslocal> >     &faceData,
		std::vector<std::vector<Element*> > &faceCodes,
		std::vector<std::vector<int> >         &faceNeighbors)
{
	size_t threads = environment->OMP_NUM_THREADS;

	for (size_t t = 1; t < threads; t++) {
		faceNeighbors[0].insert(faceNeighbors[0].end(), faceNeighbors[t].begin(), faceNeighbors[t].end());
	}
	Esutils::threadDistributionToFullDistribution(faceDistribution);
	boundaries->faces = new serializededata<eslocal, eslocal>(faceDistribution, faceData);
	boundaries->facepointers = new serializededata<eslocal, Element*>(1, faceCodes);
	std::vector<eslocal> permutation(boundaries->faces->structures());
	std::vector<size_t> fdistribution = tarray<eslocal>::distribute(threads, permutation.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		if (faceNeighbors[0][2 * i] == faceNeighbors[0][2 * j]) {
			return faceNeighbors[0][2 * i + 1] < faceNeighbors[0][2 * j + 1];
		}
		return faceNeighbors[0][2 * i] < faceNeighbors[0][2 * j];
	});
	boundaries->faces->permute(permutation, &fdistribution);

	std::vector<std::vector<EInterval> > fintervals(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = fdistribution[t]; i < fdistribution[t + 1]; ++i) {
			if (
					i > fdistribution[t] &&
					faceNeighbors[0][2 * permutation[i - 1]] == faceNeighbors[0][2 * permutation[i]] &&
					faceNeighbors[0][2 * permutation[i - 1] + 1] == faceNeighbors[0][2 * permutation[i] + 1]) {
				++fintervals[t].back().end;
			} else {
				fintervals[t].push_back(EInterval(i, i + 1, { faceNeighbors[0][2 * permutation[i]], faceNeighbors[0][2 * permutation[i] + 1] }));
			}
		}
	}

	for (size_t t = 1; t < threads; t++) {
		if (fintervals[t].size() && fintervals[0].back().neighbors.front() == fintervals[t].front().neighbors.front()) {
			fintervals[0].back().end = fintervals[t].front().end;
			fintervals[0].insert(fintervals[0].end(), fintervals[t].begin() + 1, fintervals[t].end());
		} else {
			fintervals[0].insert(fintervals[0].end(), fintervals[t].begin(), fintervals[t].end());
		}
	}
	boundaries->facesIntervals = fintervals[0];
}

void Transformation::distributeNodesToIntervals(Mesh &mesh, BoundaryStore* &boundaries)
{
	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > bnodes(boundaries->facesIntervals.size() + 1);

	#pragma omp parallel for
	for (size_t i = 0; i < boundaries->facesIntervals.size() + 1; ++i) {
		if (i == bnodes.size() - 1) {
			auto begin = boundaries->faces->datatarray().begin();
			auto end   = boundaries->faces->datatarray().end();
			bnodes[i].insert(bnodes[i].end(), begin, end);
			Esutils::sortAndRemoveDuplicity(bnodes[i]);
		} else {
			auto begin = (boundaries->faces->cbegin() + boundaries->facesIntervals[i].begin)->begin();
			auto end   = (boundaries->faces->cbegin() + boundaries->facesIntervals[i].end)->begin();
			bnodes[i].insert(bnodes[i].end(), begin, end);
			Esutils::sortAndRemoveDuplicity(bnodes[i]);
		}
	}

	std::vector<std::vector<eslocal> > nodeNeighDistribution(threads), nodeNeighData(threads);
	std::vector<size_t> ndistribution = tarray<eslocal>::distribute(threads, bnodes.back().size());

	nodeNeighDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<int> neighs;
		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			neighs.clear();
			for (size_t i = 0; i < boundaries->facesIntervals.size(); ++i) {
				if (std::binary_search(bnodes[i].begin(), bnodes[i].end(), bnodes.back()[n])) {
					neighs.insert(neighs.end(), boundaries->facesIntervals[i].neighbors.begin(), boundaries->facesIntervals[i].neighbors.end());
				}
			}
			Esutils::sortAndRemoveDuplicity(neighs);
			nodeNeighData[t].insert(nodeNeighData[t].end(), neighs.begin(), neighs.end());
			nodeNeighDistribution[t].push_back(nodeNeighData[t].size());
		}
	}

	Esutils::threadDistributionToFullDistribution(nodeNeighDistribution);
	serializededata<eslocal, eslocal> nodesNeigh(nodeNeighDistribution, nodeNeighData);

	std::vector<eslocal> permutation(nodesNeigh.structures());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		auto ni = nodesNeigh.cbegin() + i;
		auto nj = nodesNeigh.cbegin() + j;
		for (size_t n = 0; n < ni->size() && n < nj->size(); ++n) {
			if ((*ni)[n] != (*nj)[n]) {
				return (*ni)[n] < (*nj)[n];
			}
		}
		if (ni->size() != nj->size()) {
			return ni->size() < nj->size();
		}
		return bnodes.back()[i] < bnodes.back()[j];
	});

	std::vector<std::vector<eslocal> > nodesData(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = ndistribution[t]; i < ndistribution[t + 1]; ++i) {
			nodesData[t].push_back(bnodes.back()[permutation[i]]);
		}
	}

	boundaries->nodes = new serializededata<eslocal, eslocal>(1, nodesData);

	std::vector<eslocal> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (eslocal i, eslocal j) { return permutation[i] < permutation[j]; });

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = boundaries->faces->begin(t)->begin(); n != boundaries->faces->end(t)->begin(); ++n) {
			*n = backpermutation[std::lower_bound(bnodes.back().begin(), bnodes.back().end(), *n) - bnodes.back().begin()];
		}
	}

	Transformation::computeIntervals(boundaries->nodesIntervals, nodesNeigh, ndistribution, permutation);
}

void Transformation::computeProcessBoundaries(Mesh &mesh)
{
	if (mesh._processBoundaries->faces != NULL && mesh._processBoundaries->nodes != NULL) {
		return;
	}

	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::computation of process boundaries started.";

	if (mesh._nodes->elems == NULL) {
		Transformation::addLinkFromTo(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
	}

	if (mesh._elems->dual == NULL) {
		Transformation::computeDual(mesh);
	}

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<esglobal> IDBoundaries = mesh._elems->gatherElementDistrubution();

	std::vector<std::vector<eslocal> > faceDistribution(threads), faceData(threads);
	std::vector<std::vector<Element*> > faceCodes(threads);
	std::vector<std::vector<int> > faceNeighbors(threads);

	Transformation::computeBoundaries(mesh, mesh._elems->dual, 0, IDBoundaries, NULL, &faceDistribution, &faceData, &faceCodes, &faceNeighbors);

	// Transformation::distributeElementsToIntervals(mesh, mesh._processBoundaries, IDBoundaries, elementData);
	Transformation::distributeFacesToIntervals(mesh, mesh._processBoundaries, faceDistribution, faceData, faceCodes, faceNeighbors);
	Transformation::distributeNodesToIntervals(mesh, mesh._processBoundaries);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of process boundaries finished.";
}

void Transformation::computeDomainsBoundaries(Mesh &mesh)
{
	if (mesh._domainsBoundaries->faces != NULL && mesh._domainsBoundaries->nodes != NULL) {
		return;
	}

	if (mesh._domains == NULL) {
		ESINFO(TVERBOSITY) << std::string(2 * (level + 1), ' ') << "MESH::computation of domains boundaries skipped. There are no domains.";
		return;
	}

	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::computation of domains boundaries started.";

	if (mesh._nodes->elems == NULL) {
		Transformation::addLinkFromTo(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
	}

	if (mesh._elems->decomposedDual == NULL) {
		Transformation::computeDecomposedDual(mesh, TFlags::SEPARATE::ETYPES);
	}

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<esglobal> IDBoundaries = mesh._domains->gatherDomainDistribution();

	std::vector<std::vector<eslocal> > faceDistribution(threads), faceData(threads);
	std::vector<std::vector<Element*> > faceCodes(threads);
	std::vector<std::vector<int> > faceNeighbors(threads);

	Transformation::computeBoundaries(mesh, mesh._elems->decomposedDual, mesh._elems->IDs->datatarray().front(), IDBoundaries, NULL, &faceDistribution, &faceData, &faceCodes, &faceNeighbors);

	Transformation::distributeFacesToIntervals(mesh, mesh._domainsBoundaries, faceDistribution, faceData, faceCodes, faceNeighbors);
	Transformation::distributeNodesToIntervals(mesh, mesh._domainsBoundaries);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of domains boundaries finished.";
}

