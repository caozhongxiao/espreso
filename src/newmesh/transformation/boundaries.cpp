
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

	if (mesh._elems->dual == NULL) {
		Transformation::computeDual(mesh);
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

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<std::vector<eslocal> > bnodes(mesh._domains->structures());
	std::vector<std::vector<esglobal> > belements(mesh._domains->structures());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		for (size_t d = mesh._domains->datatarray().distribution()[t]; d < mesh._domains->datatarray().distribution()[t + 1]; ++d) {
			MeshDomain *domain = mesh._domains->datatarray()[d];
			auto dual = mesh._elems->dual->cbegin() + domain->eoffset;
			auto epointer = mesh._elems->epointers->cbegin() + domain->eoffset;

			std::cout << "offset: " << domain->eoffset << "\n";
			for (auto e = domain->elements->begin(); e != domain->elements->end(); ++e, ++dual, epointer) {
				if (dual->size() < epointer->front()->faces->structures()) {
					std::cout << "on boundary\n";
				} else {
					std::cout << "inner\n";
				}
			}
		}
	}

//
//	std::vector<std::vector<eslocal> > bnodes(threads);
//	std::vector<std::vector<esglobal> > bfaces(threads);
//
//	std::vector<eslocal> doffsets;
//	for (size_t d = 0; d < mesh._domains->structures(); d++) {
//		doffsets.push_back(mesh._domains->datatarray()[d]->eoffset);
//	}
//	doffsets.push_back(mesh._elems->size);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		auto elems = mesh._nodes->elems->cbegin(t);
//
//		for (size_t n = mesh._nodes->distribution[t]; n < mesh._nodes->distribution[t + 1]; ++n, ++elems) {
//			for (auto e = elems->begin(); e != elems->end(); ++e) {
//				if (mesh._elems->IDs->datatarray().front() <= *e && *e <= mesh._elems->IDs->datatarray().back()) {
//
//				}
//			}
//		}
//	}
//
//	#pragma omp parallel for
//	for (size_t t= 0; t < threads; t++) {
//		auto domains = mesh._domains->datatarray();
//		for (size_t d = domains.distribution()[t]; d < domains.distribution()[t + 1]; ++d) {
//
//
//			auto element  = mesh._domains->datatarray()[d]->elements->begin();
//			auto dual     = mesh._elems->decomposedDual->cbegin() + mesh._domains->datatarray()[d]->eoffset;
//			auto epointer = mesh._elems->epointers->datatarray();
//
//			for (size_t e = domains[d]->eoffset; e < domains[d]->eoffset + domains[d]->esize; ++e, ++element, ++dual) {
//				if (dual->size() < epointer[e]->faces->structures()) {
//					std::cout << "e: " << e << "\n";
//				} // else -> inner node
//				for (auto n = dual->begin(); n != dual->end(); ++n) {
//					if ()
//				}
//			}
//		}
//	}


	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of domains boundaries finished.";
}

