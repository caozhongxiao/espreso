
#include "transformations.h"

#include "../newmesh.h"
#include "../elements/elementstore.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"
#include "../../config/ecf/environment.h"

#include <algorithm>

using namespace espreso;

void Transformation::addLinkFromTo(NewMesh &mesh, TFlags::ELEVEL from, TFlags::ELEVEL to)
{
	auto elevel = [] (TFlags::ELEVEL level) {
		switch (level) {
		case TFlags::ELEVEL::ELEMENT:
			return "ELEMENTS";
		case TFlags::ELEVEL::FACE:
			return "FACES";
		case TFlags::ELEVEL::EDGE:
			return "EDGES";
		case TFlags::ELEVEL::NODE:
			return "NODES";
		default:
			return "??";
		}
	};
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::add links from " << elevel(from) << " to " << elevel(to) << " started.";

	ElementStore *storefrom(NULL), *storeto(NULL);
	switch (from) {
	case TFlags::ELEVEL::NODE:
		storefrom = mesh._nodes;
		break;
	default:
		break;
	}

	switch (to) {
	case TFlags::ELEVEL::ELEMENT:
		storeto = mesh._elems;
		break;
	case TFlags::ELEVEL::FACE:
		storeto = mesh._faces;
		break;
	case TFlags::ELEVEL::EDGE:
		storeto = mesh._edges;
		break;
	default:
		break;
	}

	if (storefrom == NULL || storeto == NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement addLinkFromTo(mesh, " << elevel(from) << ", " << elevel(to) << ").";
	}

	size_t threads = environment->OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(mesh._neighbours.begin(), mesh._neighbours.end(), neighbour) - mesh._neighbours.begin();
	};

	// thread x neighbor x vector(from, to)
	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > sBuffer(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(mesh._neighbours.size()));
	// neighbor x (from, to)
	std::vector<std::vector<std::pair<esglobal, esglobal> > > rBuffer(mesh._neighbours.size());
	// thread x vector(from, to)
	std::vector<std::vector<std::pair<esglobal, esglobal> > > localLinks(threads);


	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nodes = storeto->nodes->cbegin(t);
		auto IDto = storeto->IDs->datatarray();

		auto IDfrom = storefrom->IDs->datatarray();

		for (size_t e = storeto->distribution[t]; e < storeto->distribution[t + 1]; ++e, ++nodes) {
			for (size_t n = 0; n < nodes->size(); ++n) {
				localLinks[t].push_back(std::make_pair(IDfrom[(*nodes)[n]], IDto[e]));

				auto ranks = storefrom->ranks->cbegin() + (*nodes)[n];
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != environment->MPIrank) {
						sBuffer[t][n2i(*rank)].push_back(localLinks[t].back());
					}
				}
			}
		}
		std::sort(localLinks[t].begin(), localLinks[t].end());
	}

	#pragma omp parallel for
	for (size_t n = 0; n < sBuffer[0].size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
		std::sort(sBuffer[0][n].begin(), sBuffer[0][n].end());
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, mesh._neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: addLinkFromTo - exchangeUnknownSize.";
	}

	std::vector<std::vector<eslocal> > linksBoundaries(threads);
	std::vector<std::vector<esglobal> > linksData(threads);

	auto compare = [] (const std::pair<esglobal, esglobal> &position, const std::pair<esglobal, esglobal> &value) { return position.first < value.first; };
	auto addLink = [&] (const std::vector<std::pair<esglobal, esglobal> > &links, const std::pair<esglobal, esglobal> &pair, size_t t) {
		auto value = std::lower_bound(links.begin(), links.end(), pair, compare);
		for(; value != links.end() && value->first == pair.first; ++value) {
			++linksBoundaries[t].back();
			linksData[t].push_back(value->second);
		}
	};

	linksBoundaries.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::pair<esglobal, esglobal> IDpair;
		auto IDfrom = storefrom->IDs->datatarray();

		for (size_t n = storefrom->distribution[t]; n < storefrom->distribution[t + 1]; ++n) {
			IDpair.first = IDfrom[n];

			linksBoundaries[t].push_back(0);
			for (size_t r = 0; r < rBuffer.size(); r++) {
				addLink(rBuffer[r], IDpair, t);
			}
			for (size_t tt = 0; tt < threads; ++tt) {
				addLink(localLinks[tt], IDpair, t);
			}

			std::sort(linksData[t].end() - linksBoundaries[t].back(), linksData[t].end());
			if (linksBoundaries[t].size() > 1) {
				linksBoundaries[t].back() += *(linksBoundaries[t].end() - 2);
			}
		}
	}

	Esutils::threadDistributionToFullDistribution(linksBoundaries);

	switch (to) {
	case TFlags::ELEVEL::ELEMENT:
		storefrom->elems = new serializededata<eslocal, eslocal>(linksBoundaries, linksData);
		break;
	case TFlags::ELEVEL::FACE:
		storefrom->faces = new serializededata<eslocal, eslocal>(linksBoundaries, linksData);
		break;
	case TFlags::ELEVEL::EDGE:
		storefrom->edges = new serializededata<eslocal, eslocal>(linksBoundaries, linksData);
		break;
	default:
		break;
	}

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::add links from " << elevel(from) << " to " << elevel(to) << " finished.";
}

struct __haloElement__ {
	esglobal id;
	int body, material, code;
};

void Transformation::exchangeHaloElements(NewMesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::exchanging halo elements started.";

	if (mesh._nodes->elems == NULL) {
		Transformation::addLinkFromTo(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
	}

	size_t threads = environment->OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(mesh._neighbours.begin(), mesh._neighbours.end(), neighbour) - mesh._neighbours.begin();
	};

	// thread x neighbors x halo elements
	std::vector<std::vector<std::vector<__haloElement__> > > sBuffer(threads, std::vector<std::vector<__haloElement__> >(mesh._neighbours.size()));
	std::vector<std::vector<__haloElement__> > rBuffer(mesh._neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto IDs = mesh._elems->IDs->datatarray();
		auto body = mesh._elems->body->datatarray();
		auto material = mesh._elems->material->datatarray();
		auto code = mesh._elems->epointers->datatarray();
		auto nodes = mesh._elems->nodes->cbegin(t);
		std::vector<int> neighbors;
		__haloElement__ haloElement;
		haloElement.id = -1;

		for (size_t e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++nodes) {
			neighbors.clear();
			for (size_t n = 0; n < nodes->size(); ++n) {
				auto ranks = mesh._nodes->ranks->cbegin() + (*nodes)[n];
				neighbors.insert(neighbors.end(), ranks->begin(), ranks->end());
			}
			std::sort(neighbors.begin(), neighbors.end());
			Esutils::removeDuplicity(neighbors);

			if (neighbors.size() > 1) {
				haloElement.id = IDs[e];
				haloElement.body = body[e];
				haloElement.material = material[e];
				haloElement.code = code[e] - mesh._eclasses[t].data();

				for (size_t n = 0; n < neighbors.size(); n++) {
					if (neighbors[n] != environment->MPIrank) {
						sBuffer[t][n2i(neighbors[n])].push_back(haloElement);
					}
				}
			}
		}
	}

	// DUAL has to be computed before re-partition -> elements IDs has to be sorted
	#pragma omp parallel for
	for (size_t n = 0; n < sBuffer[0].size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, mesh._neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange halo elements.";
	}

	std::vector<std::vector<esglobal> > hid(threads);
	std::vector<std::vector<int> > hbody(threads), hmaterial(threads);
	std::vector<std::vector<NewElement*> > hcode(threads);

	for (size_t n = 0; n < rBuffer.size(); ++n) {
		std::vector<size_t> distribution = tarray<esglobal>::distribute(threads, rBuffer[n].size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; ++e) {
				hid[t].push_back(rBuffer[n][e].id);
				hbody[t].push_back(rBuffer[n][e].body);
				hmaterial[t].push_back(rBuffer[n][e].material);
				hcode[t].push_back(mesh._eclasses[0].data() + rBuffer[n][e].code);
			}
		}
	}

	serializededata<eslocal, esglobal>::balance(1, hid);
	serializededata<eslocal, eslocal>::balance(1, hbody);
	serializededata<eslocal, eslocal>::balance(1, hmaterial);
	serializededata<eslocal, NewElement*>::balance(1, hcode);

	mesh._halo->IDs = new serializededata<eslocal, esglobal>(1, hid);
	mesh._halo->body = new serializededata<eslocal, eslocal>(1, hbody);
	mesh._halo->material = new serializededata<eslocal, eslocal>(1, hmaterial);
	mesh._halo->epointers = new serializededata<eslocal, NewElement*>(1, hcode);

	mesh._halo->size = mesh._halo->IDs->datatarray().size();
	mesh._halo->distribution = mesh._halo->IDs->datatarray().distribution();

	mesh._halo->sort();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto epointer = mesh._halo->epointers->datatarray();

		for (auto e = mesh._halo->distribution[t]; e < mesh._halo->distribution[t + 1]; ++e) {
			epointer[e] = mesh._eclasses[t].data() + (epointer[e] - mesh._eclasses[0].data());
		}
	}

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::exchanging halo elements finished.";
}

void Transformation::computeDual(NewMesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::computation of the dual graph of all elements started.";

	if (mesh._nodes->elems == NULL) {
		Transformation::addLinkFromTo(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
	}

	if (mesh._halo->IDs == NULL) {
		Transformation::exchangeHaloElements(mesh);
	}

	size_t threads = environment->OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(mesh._neighbours.begin(), mesh._neighbours.end(), neighbour) - mesh._neighbours.begin();
	};

	std::vector<std::vector<eslocal> > dualDistribution(threads);
	std::vector<std::vector<esglobal> > dualData(threads);

	dualDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto IDs = mesh._elems->IDs->datatarray();
		auto nodes = mesh._elems->nodes->cbegin(t);
		std::vector<esglobal> neighElementIDs;
		int myCommon, neighCommon;
		size_t neigh, nCounter;

		for (size_t e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++nodes) {
			neighElementIDs.clear();
			for (size_t n = 0; n < nodes->size(); ++n) {
				auto elements = mesh._nodes->elems->cbegin() + (*nodes)[n];
				neighElementIDs.insert(neighElementIDs.end(), elements->begin(), elements->end());
			}
			std::sort(neighElementIDs.begin(), neighElementIDs.end());
			myCommon = mesh._elems->epointers->datatarray()[e]->nCommonFace;

			neigh = 0;
			nCounter = 1;
			dualDistribution[t].push_back(0);
			while (neigh < neighElementIDs.size()) {
				while (neigh + nCounter < neighElementIDs.size() && neighElementIDs[neigh] == neighElementIDs[neigh + nCounter]) {
					nCounter++;
				}

				if (IDs[e] != neighElementIDs[neigh]) {
					auto it = std::lower_bound(IDs.begin(), IDs.end(), neighElementIDs[neigh]);
					if (it != IDs.end()) {
						neighCommon = mesh._elems->epointers->datatarray()[it - IDs.begin()]->nCommonFace;
					} else {
						auto it = std::lower_bound(mesh._halo->IDs->datatarray().begin(), mesh._halo->IDs->datatarray().end(), neighElementIDs[neigh]);
						neighCommon = mesh._halo->epointers->datatarray()[it - mesh._halo->IDs->datatarray().begin()]->nCommonFace;
					}
					if (nCounter >= std::min(myCommon, neighCommon)) {
						++dualDistribution[t].back();
						dualData[t].push_back(neighElementIDs[neigh]);
					}
				}
				neigh += nCounter;
				nCounter = 1;
			}
			if (dualDistribution[t].size() > 1) {
				dualDistribution[t].back() += *(dualDistribution[t].end() - 2);
			}
		}
	}

	Esutils::threadDistributionToFullDistribution(dualDistribution);

	mesh._elems->dual = new serializededata<eslocal, esglobal>(dualDistribution, dualData);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of the dual graph of all elements finished.";
}

void Transformation::computeDecomposedDual(NewMesh &mesh, TFlags::SEPARATE separate)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::computation of the decomposed dual graph of local elements started.";

	if (mesh._nodes->elems == NULL) {
		Transformation::addLinkFromTo(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > dualDistribution(threads);
	std::vector<std::vector<esglobal> > dualData(threads);

	dualDistribution.front().push_back(0);

	if (mesh._elems->dual == NULL) {
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto IDs = mesh._elems->IDs->datatarray();
			auto nodes = mesh._elems->nodes->cbegin(t);
			auto body = mesh._elems->body->datatarray();
			auto material = mesh._elems->body->datatarray();
			auto epointer = mesh._elems->epointers->datatarray();
			std::vector<esglobal> neighElementIDs;
			int myCommon, neighCommon;
			size_t neigh, nCounter;

			for (size_t e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++nodes) {
				neighElementIDs.clear();
				for (size_t n = 0; n < nodes->size(); ++n) {
					auto elements = mesh._nodes->elems->cbegin() + (*nodes)[n];
					neighElementIDs.insert(neighElementIDs.end(), elements->begin(), elements->end());
				}
				std::sort(neighElementIDs.begin(), neighElementIDs.end());
				myCommon = mesh._elems->epointers->datatarray()[e]->nCommonFace;

				neigh = 0;
				nCounter = 1;
				dualDistribution[t].push_back(0);
				while (neigh < neighElementIDs.size()) {
					while (neigh + nCounter < neighElementIDs.size() && neighElementIDs[neigh] == neighElementIDs[neigh + nCounter]) {
						nCounter++;
					}

					if (IDs[e] != neighElementIDs[neigh]) {
						auto it = std::lower_bound(IDs.begin(), IDs.end(), neighElementIDs[neigh]);
						if (it != IDs.end()) {
							neighCommon = mesh._elems->epointers->datatarray()[it - IDs.begin()]->nCommonFace;
							if (
									nCounter >= std::min(myCommon, neighCommon) &&
									body[e] == body[it - IDs.begin()] &&
									(!(separate & TFlags::SEPARATE::MATERIALS) || material[e] == material[it - IDs.begin()]) &&
									(!(separate & TFlags::SEPARATE::ETYPES) || epointer[e]->type == epointer[it - IDs.begin()]->type)
								) {


								if (
										(e  < mesh._elems->size / 2 && it - IDs.begin()  < mesh._elems->size / 2) ||
										(e >= mesh._elems->size / 2 && it - IDs.begin() >= mesh._elems->size / 2)) {

									++dualDistribution[t].back();
									dualData[t].push_back(neighElementIDs[neigh]);
								}
							}
						}
					}
					neigh += nCounter;
					nCounter = 1;
				}
				if (dualDistribution[t].size() > 1) {
					dualDistribution[t].back() += *(dualDistribution[t].end() - 2);
				}
			}
		}
	} else {
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto IDs = mesh._elems->IDs->datatarray();
			auto dual = mesh._elems->dual->cbegin(t);
			auto body = mesh._elems->body->datatarray();
			auto material = mesh._elems->body->datatarray();
			auto epointer = mesh._elems->epointers->datatarray();

			for (size_t e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++dual) {
				dualDistribution[t].push_back(0);
				for (auto neigh = dual->begin(); neigh != dual->end(); ++neigh) {
					auto it = std::lower_bound(IDs.begin(), IDs.end(), *neigh);
					if (
							it != IDs.end() &&
							body[e] == body[it - IDs.begin()] &&
							(!(separate & TFlags::SEPARATE::MATERIALS) || material[e] == material[it - IDs.begin()]) &&
							(!(separate & TFlags::SEPARATE::ETYPES) || epointer[e]->type == epointer[it - IDs.begin()]->type)
						) {

						++dualDistribution[t].back();
						dualData[t].push_back(*neigh);
					}
				}
				if (dualDistribution[t].size() > 1) {
					dualDistribution[t].back() += *(dualDistribution[t].end() - 2);
				}
			}
		}
	}

	Esutils::threadDistributionToFullDistribution(dualDistribution);

	mesh._elems->decomposedDual = new serializededata<eslocal, esglobal>(dualDistribution, dualData);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::computation of the decomposed dual graph of local elements finished.";
}



