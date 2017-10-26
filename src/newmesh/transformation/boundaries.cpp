
#include "transformations.h"

#include "../newmesh.h"
#include "../elements/newelement.h"
#include "../elements/elementstore.h"
#include "../store/domainstore.h"
#include "../store/boundarystore.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../config/ecf/environment.h"

#include <algorithm>
#include <iostream>

using namespace espreso;

void Transformation::computeProcessBoundaries(NewMesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::computation of process boundaries started.";

	if (mesh._nodes->elems == NULL) {
		Transformation::addLinkFromTo(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
	}

	if (mesh._elems->dual == NULL) {
		Transformation::computeDual(mesh);
	}

	std::vector<esglobal> IDBoundaries = mesh._elems->gatherSizes();
	esglobal begine = IDBoundaries[environment->MPIrank];
	esglobal ende   = IDBoundaries[environment->MPIrank + 1];

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > elementData(threads), faceDistribution(threads), faceData(threads), faceNodes(threads);
	std::vector<std::vector<NewElement*> > faceCodes(threads);
	std::vector<std::vector<int> > neighProcess(threads);

	faceDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esglobal> common;
		size_t ncommons, counter, neighbor;
		auto dual = mesh._elems->dual->cbegin(t);
		auto epointer = mesh._elems->epointers->cbegin(t);
		esglobal eID = mesh._elems->distribution[t];

		for (auto e = mesh._elems->nodes->cbegin(t); e != mesh._elems->nodes->cend(t); ++e, ++dual, ++epointer, ++eID) {
			if (dual->front() < begine || dual->back() >= ende || dual->size() < epointer->front()->faces->structures()) {

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
						faceDistribution[t].push_back(face->size());
						faceCodes[t].push_back(facepointer->front());
						if (faceDistribution[t].size() > 1) {
							faceDistribution[t].back() += *(faceDistribution[t].end() - 2);
						}
						for (auto n = face->begin(); n != face->end(); ++n) {
							faceData[t].push_back((*e)[*n]);
						}
						neighProcess[t].push_back(neighbor);
						if (neighbor != -1) {
							elementData[t].push_back(eID);
						}
					}
				}
			}
		}
		faceNodes[t] = faceData[t];
		Esutils::sortAndRemoveDuplicity(faceNodes[t]);
		Esutils::sortAndRemoveDuplicity(elementData[t]);
	}
	Esutils::mergeThreadedUniqueData(faceNodes);
	for (size_t t = 1; t < threads; t++) {
		neighProcess[0].insert(neighProcess[0].end(), neighProcess[t].begin(), neighProcess[t].end());
	}
	faceNodes.resize(1);
	faceNodes.resize(threads);

	Esutils::threadDistributionToFullDistribution(faceDistribution);

	serializededata<eslocal, eslocal>::balance(1, elementData);
	mesh._processBoundaries->elems = new serializededata<eslocal, eslocal>(1, elementData);

	mesh._processBoundaries->faces = new serializededata<eslocal, eslocal>(faceDistribution, faceData);
	mesh._processBoundaries->facepointers = new serializededata<eslocal, NewElement*>(1, faceCodes);
	std::vector<eslocal> permutation(mesh._processBoundaries->faces->structures());
	std::vector<size_t> fdistribution = tarray<eslocal>::distribute(threads, permutation.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return neighProcess[i] < neighProcess[j]; });
	mesh._processBoundaries->faces->permute(permutation, &fdistribution);

	mesh._processBoundaries->nodes = new serializededata<eslocal, eslocal>(1, faceNodes);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = mesh._processBoundaries->faces->begin(t)->begin(); n != mesh._processBoundaries->faces->end(t)->begin(); ++n) {
			*n = std::lower_bound(faceNodes[0].begin(), faceNodes[0].end(), *n) - faceNodes[0].begin();
		}
	}

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of process boundaries finished.";
}

void Transformation::computeDomainsBoundaries(NewMesh &mesh)
{
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

	std::vector<esglobal> IDBoundaries = mesh._elems->gatherSizes();

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<std::vector<eslocal> > faceDistribution(threads), faceData(threads), faceNodes(threads);
	std::vector<std::vector<NewElement*> > faceCodes(threads);

	faceDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esglobal> common;
		size_t ncommons, counter;
		bool isAdept;

		for (size_t d = mesh._domains->domainDistribution[t]; d < mesh._domains->domainDistribution[t + 1]; ++d) {
			esglobal begine = IDBoundaries[environment->MPIrank] + mesh._domains->domainElementBoundaries[d];
			esglobal ende   = IDBoundaries[environment->MPIrank] + mesh._domains->domainElementBoundaries[d + 1];
			auto dual = mesh._elems->decomposedDual->cbegin() + mesh._domains->domainElementBoundaries[d];
			auto epointer = mesh._elems->epointers->cbegin() + mesh._domains->domainElementBoundaries[d];

			for (
					auto e = mesh._elems->nodes->cbegin() + mesh._domains->domainElementBoundaries[d];
					e != mesh._elems->nodes->cbegin() + mesh._domains->domainElementBoundaries[d + 1];
					++e, ++dual, ++epointer) {

				isAdept = false;
				if (dual->size() < epointer->front()->faces->structures()) {
					isAdept = true;
				} else {
					for (auto ne = dual->begin(); ne != dual->end(); ++ne) {
						if (*ne < begine || ende <= *ne) {
							isAdept = true;
							break;
						}
					}
				}

				if (isAdept) {
					auto facepointer = epointer->front()->facepointers->cbegin(t);
					for (auto face = epointer->front()->faces->cbegin(t); face != epointer->front()->faces->cend(t); ++face, ++facepointer) {

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
									}
								}
								counter = 0;
							}
						}
						if (face->size() == counter + 1) {
							if (begine <= common.back() && common.back() < ende) {
								++ncommons;
							}
						}

						if (ncommons == 1) {
							faceDistribution[t].push_back(face->size());
							faceCodes[t].push_back(facepointer->front());
							if (faceDistribution[t].size() > 1) {
								faceDistribution[t].back() += *(faceDistribution[t].end() - 2);
							}
							for (auto n = face->begin(); n != face->end(); ++n) {
								faceData[t].push_back((*e)[*n]);
							}
						}
					}
				}
			}
		}
		faceNodes[t] = faceData[t];
		Esutils::sortAndRemoveDuplicity(faceNodes[t]);
	}
	Esutils::mergeThreadedUniqueData(faceNodes);
	faceNodes.resize(1);
	faceNodes.resize(threads);

	Esutils::threadDistributionToFullDistribution(faceDistribution);

	mesh._domainsBoundaries->faces = new serializededata<eslocal, eslocal>(faceDistribution, faceData);
	mesh._domainsBoundaries->facepointers = new serializededata<eslocal, NewElement*>(1, faceCodes);
	mesh._domainsBoundaries->nodes = new serializededata<eslocal, eslocal>(1, faceNodes);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = mesh._domainsBoundaries->faces->begin(t)->begin(); n != mesh._domainsBoundaries->faces->end(t)->begin(); ++n) {
			*n = std::lower_bound(faceNodes[0].begin(), faceNodes[0].end(), *n) - faceNodes[0].begin();
		}
	}

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of domains boundaries finished.";
}

