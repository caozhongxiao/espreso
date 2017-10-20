
#include "transformations.h"

#include "../newmesh.h"
#include "../elements/newelement.h"
#include "../elements/elementstore.h"
#include "../store/boundarystore.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../config/ecf/environment.h"

#include <algorithm>
#include <iostream>

using namespace espreso;

void Transformation::computeProcessesCommonBoundary(NewMesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::computation of processes boundary started.";

	if (mesh._nodes->elems == NULL) {
		Transformation::addLinkFromTo(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
	}

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<std::vector<eslocal> > bnodes(threads);
	std::vector<std::vector<esglobal> > belements(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = mesh._nodes->ranks->cbegin(t);
		auto elems = mesh._nodes->elems->cbegin(t);

		for (size_t n = mesh._nodes->distribution[t]; n < mesh._nodes->distribution[t + 1]; ++n, ++ranks, ++elems) {
			if (ranks->size() > 1) {
				bnodes[t].push_back(n);
				for (auto e = elems->begin(); e != elems->end(); ++e) {
					if (mesh._elems->IDs->datatarray().front() <= *e && *e <= mesh._elems->IDs->datatarray().back()) {
						belements[t].push_back(*e - mesh._elems->IDs->datatarray().front());
					}
				}
			}
		}
		Esutils::sortAndRemoveDuplicity(belements[t]);
	}
	Esutils::mergeThreadedUniqueData(belements);
	belements.resize(1);
	belements.resize(threads);

	serializededata<eslocal, esglobal>::balance(1, belements);
	serializededata<eslocal, eslocal>::balance(1, bnodes);

	mesh._processesCommonBoundary->elems = new serializededata<eslocal, esglobal>(1, belements);
	mesh._processesCommonBoundary->nodes = new serializededata<eslocal, eslocal>(1, bnodes);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of processes boundary finished.";
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

	if (mesh._nodes->domains == NULL) {
		Transformation::assignDomainsToNodes(mesh);
	}

//	if (elevel & TFlags::ELEVEL::ELEMENT) {
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement domains boundaries for elements.";
//	}
//
//	if (elevel & TFlags::ELEVEL::EDGE) {
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement domains boundaries for elements.";
//	}

	std::vector<esglobal> IDBoundaries = mesh._elems->gatherSizes();

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<std::vector<eslocal> > faceDistribution(threads), faceData(threads), faceNodes(threads);
	std::vector<std::vector<NewElement*> > faceCodes(threads);

	faceDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto elems = mesh._nodes->elems->cbegin();
		std::vector<esglobal> common;
		size_t ncommons, counter;
		bool isAdept;

		for (size_t d = mesh._domains->datatarray().distribution()[t]; d < mesh._domains->datatarray().distribution()[t + 1]; ++d) {
			MeshDomain *domain = mesh._domains->datatarray()[d];
			esglobal begine = IDBoundaries[environment->MPIrank] + domain->eoffset;
			esglobal ende   = IDBoundaries[environment->MPIrank] + domain->eoffset + domain->esize;
			auto dual = mesh._elems->decomposedDual->cbegin() + domain->eoffset;
			auto epointer = mesh._elems->epointers->cbegin() + domain->eoffset;

			for (auto e = domain->elements->begin(); e != domain->elements->end(); ++e, ++dual, ++epointer) {

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
					for (auto face = epointer->front()->faces->cbegin(); face != epointer->front()->faces->cend(); ++face, ++facepointer) {

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

	mesh._domainsBoundaries = new BoundaryStore();
	mesh._domainsBoundaries->clusterfaces = new serializededata<eslocal, eslocal>(faceDistribution, faceData);
	mesh._domainsBoundaries->localfaces = new serializededata<eslocal, eslocal>(faceDistribution, faceData);
	mesh._domainsBoundaries->facepointers = new serializededata<eslocal, NewElement*>(1, faceCodes);
	mesh._domainsBoundaries->nodes = new serializededata<eslocal, eslocal>(1, faceNodes);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = mesh._domainsBoundaries->localfaces->begin(t)->begin(); n != mesh._domainsBoundaries->localfaces->end(t)->begin(); ++n) {
			*n = std::lower_bound(faceNodes[0].begin(), faceNodes[0].end(), *n) - faceNodes[0].begin();
		}
	}

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of domains boundaries finished.";
}

