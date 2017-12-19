
#include "meshpreprocessing.h"

#include "../mesh.h"

#include "../store/store.h"
#include "../store/elementstore.h"
#include "../store/nodestore.h"
#include "../store/elementsregionstore.h"
#include "../store/boundaryregionstore.h"
#include "../store/sharedinterfacestore.h"

#include "../elements/element.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/parser.h"
#include "../../basis/logging/logging.h"

#include "../../config/ecf/environment.h"

#include "../../wrappers/wparmetis.h"
#include "../../wrappers/wmetis.h"

#include <algorithm>
#include <numeric>
#include <cstring>

using namespace espreso;

size_t MeshPreprocessing::level = 0;

void MeshPreprocessing::start(const std::string &message)
{
	ESINFO(VERBOSITY(level)) << std::string(2 * level, ' ') << "Mesh preprocessing :: " << message << " started.";
	++level;
}

void MeshPreprocessing::skip(const std::string &message)
{
	ESINFO(VERBOSITY(level)) << std::string(2 * level, ' ') << "Mesh preprocessing :: " << message << " skipped.";
}

void MeshPreprocessing::finish(const std::string &message)
{
	--level;
	ESINFO(VERBOSITY(level)) << std::string(2 * level, ' ') << "Mesh preprocessing :: " << message << " finished.";
}

void MeshPreprocessing::linkNodesAndElements()
{
	start("link nodes and elements");

	if (_mesh->elements == NULL || _mesh->nodes == NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: fill both elements and nodes.";
	}

	size_t threads = environment->OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbour) - _mesh->neighbours.begin();
	};

	// thread x neighbor x vector(from, to)
	std::vector<std::vector<std::vector<std::pair<eslocal, eslocal> > > > sBuffer(threads, std::vector<std::vector<std::pair<eslocal, eslocal> > >(_mesh->neighbours.size()));
	// neighbor x (from, to)
	std::vector<std::vector<std::pair<eslocal, eslocal> > > rBuffer(_mesh->neighbours.size());
	// thread x vector(from, to)
	std::vector<std::vector<std::pair<eslocal, eslocal> > > localLinks(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nodes = _mesh->elements->nodes->cbegin(t);
		const auto &IDto = _mesh->elements->IDs->datatarray();

		const auto &IDfrom = _mesh->nodes->IDs->datatarray();

		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++nodes) {
			for (size_t n = 0; n < nodes->size(); ++n) {
				localLinks[t].push_back(std::make_pair(IDfrom[(*nodes)[n]], IDto[e]));

				auto ranks = _mesh->nodes->ranks->cbegin() + (*nodes)[n];
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

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, _mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: addLinkFromTo - exchangeUnknownSize.";
	}

	std::vector<std::vector<eslocal> > linksBoundaries(threads);
	std::vector<std::vector<eslocal> > linksData(threads);

	auto compare = [] (const std::pair<eslocal, eslocal> &position, const std::pair<eslocal, eslocal> &value) { return position.first < value.first; };
	auto addLink = [&] (const std::vector<std::pair<eslocal, eslocal> > &links, const std::pair<eslocal, eslocal> &pair, size_t t) {
		auto value = std::lower_bound(links.begin(), links.end(), pair, compare);
		for(; value != links.end() && value->first == pair.first; ++value) {
			++linksBoundaries[t].back();
			linksData[t].push_back(value->second);
		}
	};

	linksBoundaries.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::pair<eslocal, eslocal> IDpair;
		const auto &IDfrom = _mesh->nodes->IDs->datatarray();

		for (size_t n = _mesh->nodes->distribution[t]; n < _mesh->nodes->distribution[t + 1]; ++n) {
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

	_mesh->nodes->elements = new serializededata<eslocal, eslocal>(linksBoundaries, linksData);

	finish("link nodes and elements");
}

struct __haloElement__ {
	eslocal id;
	int body, material, code;
};

void MeshPreprocessing::exchangeHalo()
{
	start("exchanging halo");

	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	size_t threads = environment->OMP_NUM_THREADS;
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbour) - _mesh->neighbours.begin();
	};

	// thread x neighbors x halo elements
	std::vector<std::vector<std::vector<__haloElement__> > > sBuffer(threads, std::vector<std::vector<__haloElement__> >(_mesh->neighbours.size()));
	std::vector<std::vector<__haloElement__> > rBuffer(_mesh->neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = _mesh->elements->IDs->datatarray();
		const auto &body = _mesh->elements->body->datatarray();
		const auto &material = _mesh->elements->material->datatarray();
		const auto &code = _mesh->elements->epointers->datatarray();
		auto nodes =_mesh->elements->nodes->cbegin(t);
		std::vector<int> neighbors;
		__haloElement__ haloElement;
		haloElement.id = -1;

		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++nodes) {
			neighbors.clear();
			for (size_t n = 0; n < nodes->size(); ++n) {
				auto ranks =_mesh->nodes->ranks->cbegin() + (*nodes)[n];
				neighbors.insert(neighbors.end(), ranks->begin(), ranks->end());
			}
			std::sort(neighbors.begin(), neighbors.end());
			Esutils::removeDuplicity(neighbors);

			if (neighbors.size() > 1) {
				haloElement.id = IDs[e];
				haloElement.body = body[e];
				haloElement.material = material[e];
				haloElement.code = code[e] - _mesh->_eclasses[t];

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

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, _mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange halo elements.";
	}

	std::vector<std::vector<eslocal> > hid(threads);
	std::vector<std::vector<int> > hbody(threads), hmaterial(threads);
	std::vector<std::vector<Element*> > hcode(threads);

	for (size_t n = 0; n < rBuffer.size(); ++n) {
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, rBuffer[n].size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; ++e) {
				hid[t].push_back(rBuffer[n][e].id);
				hbody[t].push_back(rBuffer[n][e].body);
				hmaterial[t].push_back(rBuffer[n][e].material);
				hcode[t].push_back(_mesh->_eclasses[0] + rBuffer[n][e].code);
			}
		}
	}

	serializededata<eslocal, eslocal>::balance(1, hid);
	serializededata<eslocal, eslocal>::balance(1, hbody);
	serializededata<eslocal, eslocal>::balance(1, hmaterial);
	serializededata<eslocal, Element*>::balance(1, hcode);

	_mesh->halo->IDs = new serializededata<eslocal, eslocal>(1, hid);
	_mesh->halo->body = new serializededata<eslocal, eslocal>(1, hbody);
	_mesh->halo->material = new serializededata<eslocal, eslocal>(1, hmaterial);
	_mesh->halo->epointers = new serializededata<eslocal, Element*>(1, hcode);

	_mesh->halo->size = _mesh->halo->IDs->datatarray().size();
	_mesh->halo->distribution = _mesh->halo->IDs->datatarray().distribution();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto &epointer = _mesh->halo->epointers->datatarray();
		for (auto e = _mesh->halo->distribution[t]; e < _mesh->halo->distribution[t + 1]; ++e) {
			epointer[e] = _mesh->_eclasses[t] + (epointer[e] - _mesh->_eclasses[0]);
		}
	}

	const auto &hIDs = _mesh->halo->IDs->datatarray();
	std::vector<eslocal> permutation(hIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return hIDs[i] < hIDs[j]; });
	_mesh->halo->permute(permutation);

	finish("exchanging halo");
}


void MeshPreprocessing::computeDual()
{
	start("computation of the dual graph of all elements");

	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	if (_mesh->halo->IDs == NULL) {
		this->exchangeHalo();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<std::vector<eslocal> > dualDistribution(threads);
	std::vector<std::vector<eslocal> > dualData(threads);

	dualDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = _mesh->elements->IDs->datatarray();
		auto nodes = _mesh->elements->nodes->cbegin(t);
		std::vector<eslocal> neighElementIDs;
		int myCommon, neighCommon;
		int neigh, nCounter;

		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++nodes) {
			neighElementIDs.clear();
			for (size_t n = 0; n < nodes->size(); ++n) {
				auto elements = _mesh->nodes->elements->cbegin() + (*nodes)[n];
				neighElementIDs.insert(neighElementIDs.end(), elements->begin(), elements->end());
			}
			std::sort(neighElementIDs.begin(), neighElementIDs.end());
			myCommon = _mesh->elements->epointers->datatarray()[e]->nCommonFace;

			neigh = 0;
			nCounter = 1;
			dualDistribution[t].push_back(0);
			while (neigh < (int)neighElementIDs.size()) {
				while (neigh + nCounter < (int)neighElementIDs.size() && neighElementIDs[neigh] == neighElementIDs[neigh + nCounter]) {
					nCounter++;
				}

				if (IDs[e] != neighElementIDs[neigh]) {
					auto it = std::lower_bound(IDs.begin(), IDs.end(), neighElementIDs[neigh]);
					if (it != IDs.end() && *it == neighElementIDs[neigh]) {
						neighCommon = _mesh->elements->epointers->datatarray()[it - IDs.begin()]->nCommonFace;
					} else {
						auto it = std::lower_bound(_mesh->halo->IDs->datatarray().begin(), _mesh->halo->IDs->datatarray().end(), neighElementIDs[neigh]);
						neighCommon = _mesh->halo->epointers->datatarray()[it - _mesh->halo->IDs->datatarray().begin()]->nCommonFace;
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

	_mesh->elements->dual = new serializededata<eslocal, eslocal>(dualDistribution, dualData);

	finish("computation of the dual graph of all elements");
}

void MeshPreprocessing::computeDecomposedDual(bool separateMaterials, bool separateEtype)
{
	start("computation of the decomposed dual graph of local elements");

	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<eslocal> IDBoundaries = _mesh->elements->gatherElementsProcDistribution();

	std::vector<std::vector<eslocal> > dualDistribution(threads);
	std::vector<std::vector<eslocal> > dualData(threads);

	dualDistribution.front().push_back(0);

	if (_mesh->elements->dual == NULL) {
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			const auto &IDs = _mesh->elements->IDs->datatarray();
			const auto &body = _mesh->elements->body->datatarray();
			const auto &material = _mesh->elements->body->datatarray();
			const auto &epointer = _mesh->elements->epointers->datatarray();
			auto nodes = _mesh->elements->nodes->cbegin(t);
			std::vector<eslocal> neighElementIDs;
			int myCommon, neighCommon;
			int neigh, nCounter;

			for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++nodes) {
				neighElementIDs.clear();
				for (size_t n = 0; n < nodes->size(); ++n) {
					auto elements = _mesh->nodes->elements->cbegin() + (*nodes)[n];
					neighElementIDs.insert(neighElementIDs.end(), elements->begin(), elements->end());
				}
				std::sort(neighElementIDs.begin(), neighElementIDs.end());
				myCommon = _mesh->elements->epointers->datatarray()[e]->nCommonFace;

				neigh = 0;
				nCounter = 1;
				dualDistribution[t].push_back(0);
				while (neigh < (eslocal)neighElementIDs.size()) {
					while (neigh + nCounter < (eslocal)neighElementIDs.size() && neighElementIDs[neigh] == neighElementIDs[neigh + nCounter]) {
						nCounter++;
					}

					if (IDs[e] != neighElementIDs[neigh]) {
						auto it = std::lower_bound(IDs.begin(), IDs.end(), neighElementIDs[neigh]);
						if (it != IDs.end() && *it == neighElementIDs[neigh]) {
							neighCommon = _mesh->elements->epointers->datatarray()[it - IDs.begin()]->nCommonFace;
							if (
									nCounter >= std::min(myCommon, neighCommon) &&
									body[e] == body[it - IDs.begin()] &&
									(!separateMaterials || material[e] == material[it - IDs.begin()]) &&
									(!separateEtype || epointer[e]->type == epointer[it - IDs.begin()]->type)
								) {

								++dualDistribution[t].back();
								dualData[t].push_back(neighElementIDs[neigh] - IDBoundaries[environment->MPIrank]);
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
			const auto &IDs = _mesh->elements->IDs->datatarray();
			const auto &body = _mesh->elements->body->datatarray();
			const auto &material = _mesh->elements->body->datatarray();
			const auto &epointer = _mesh->elements->epointers->datatarray();
			auto dual = _mesh->elements->dual->cbegin(t);

			for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++dual) {
				dualDistribution[t].push_back(0);
				for (auto neigh = dual->begin(); neigh != dual->end(); ++neigh) {
					auto it = std::lower_bound(IDs.begin(), IDs.end(), *neigh);
					if (
							it != IDs.end() && *it == *neigh &&
							body[e] == body[it - IDs.begin()] &&
							(!separateMaterials || material[e] == material[it - IDs.begin()]) &&
							(!separateEtype || epointer[e]->type == epointer[it - IDs.begin()]->type)
						) {

						++dualDistribution[t].back();
						dualData[t].push_back(*neigh - IDBoundaries[environment->MPIrank]);
					}
				}
				if (dualDistribution[t].size() > 1) {
					dualDistribution[t].back() += *(dualDistribution[t].end() - 2);
				}
			}
		}
	}

	Esutils::threadDistributionToFullDistribution(dualDistribution);

	_mesh->elements->decomposedDual = new serializededata<eslocal, eslocal>(dualDistribution, dualData);

	finish("computation of the decomposed dual graph of local elements");
}

void MeshPreprocessing::computeBoundaryNodes(std::vector<eslocal> &externalBoundary, std::vector<eslocal> &internalBoundary)
{
	start("computation of boundary nodes");

	if (_mesh->elements->dual == NULL) {
		this->computeDual();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<eslocal> IDBoundaries = _mesh->elements->gatherElementsDistribution();

	std::vector<std::vector<eslocal> > external(threads), internal(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> common;
		size_t ncommons, counter;
		bool isExternal;
		eslocal eID = _mesh->elements->distribution[t], eoffset = _mesh->elements->IDs->datatarray().front();
		auto dual = _mesh->elements->dual->cbegin(t);
		auto epointer = _mesh->elements->epointers->cbegin(t);
		auto IDpointer = std::lower_bound(IDBoundaries.begin(), IDBoundaries.end(), eID + eoffset + 1) - 1;
		eslocal begine = *IDpointer, ende = *(IDpointer + 1);

		for (auto e = _mesh->elements->nodes->cbegin(t); e != _mesh->elements->nodes->cend(t); ++e, ++dual, ++epointer, ++eID) {
			if (eID + eoffset >= ende) {
				++IDpointer;
				begine = *IDpointer;
				ende = *(IDpointer + 1);
			}
			if (dual->size() < epointer->front()->faces->structures() || dual->front() < begine || dual->back() >= ende) {

				auto facepointer = epointer->front()->facepointers->cbegin();
				for (auto face = epointer->front()->faces->cbegin(); face != epointer->front()->faces->cend(); ++face, ++facepointer) {

					isExternal = true;
					common.clear();
					for (auto n = face->begin(); n != face->end(); ++n) {
						auto nelements = _mesh->nodes->elements->cbegin() + (*e)[*n];
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
									isExternal = false;
								}
							}
							counter = 0;
						}
					}
					if (face->size() == counter + 1) {
						if (begine <= common.back() && common.back() < ende) {
							++ncommons;
						} else {
							isExternal = false;
						}
					}

					if (ncommons == 1) {
						if (isExternal) {
							for (auto n = face->begin(); n != face->end(); ++n) {
								external[t].push_back((*e)[*n]);
							}
						} else {
							for (auto n = face->begin(); n != face->end(); ++n) {
								internal[t].push_back((*e)[*n]);
							}
						}
					}
				}
			}
		}
		Esutils::sortAndRemoveDuplicity(internal[t]);
		Esutils::sortAndRemoveDuplicity(external[t]);
	}

	for (size_t t = 0; t < threads; t++) {
		externalBoundary.insert(externalBoundary.end(), external[t].begin(), external[t].end());
	}
	Esutils::sortAndRemoveDuplicity(externalBoundary);

	for (size_t t = 1; t < threads; t++) {
		internal[0].insert(internal[0].end(), internal[t].begin(), internal[t].end());
	}
	Esutils::sortAndRemoveDuplicity(internal[0]);

	internalBoundary.resize(internal[0].size());
	internalBoundary.resize(std::set_difference(internal[0].begin(), internal[0].end(), externalBoundary.begin(), externalBoundary.end(), internalBoundary.begin()) - internalBoundary.begin());

	finish("computation of boundary nodes");
}

void MeshPreprocessing::computeRegionArea(BoundaryRegionStore *store)
{
	double A = 0;
	auto nodes = store->elements->cbegin();
	const auto &epointers = store->epointers->datatarray();
	const auto &coordinates = _mesh->nodes->coordinates->datatarray();
	for (size_t e = 0; e < store->elements->structures(); ++e, ++nodes) {

		DenseMatrix coords(nodes->size(), 3), dND(1, 3);

		const std::vector<DenseMatrix> &dN = *epointers[e]->dN;
		const std::vector<double> &weighFactor = *epointers[e]->weighFactor;

		for (size_t n = 0; n < nodes->size(); ++n) {
			coords(n, 0) = coordinates[nodes->at(n)].x;
			coords(n, 1) = coordinates[nodes->at(n)].y;
			coords(n, 2) = coordinates[nodes->at(n)].z;
		}

		if (store->dimension == 1) {
			for (size_t gp = 0; gp < dN.size(); gp++) {
				dND.multiply(dN[gp], coords);
				A += dND.norm() * weighFactor[gp];
			}
		}
		if (store->dimension == 2) {
			for (size_t gp = 0; gp < dN.size(); gp++) {
				dND.multiply(dN[gp], coords);
				Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
				Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
				Point va = Point::cross(v1, v2);
				A += va.norm() * weighFactor[gp];
			}
		}
	}

	MPI_Allreduce(&A, &store->area, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
}

void MeshPreprocessing::computeSharedFaces()
{
	start("computation of shared faces");

	if (_mesh->nodes->elements == NULL) {
		linkNodesAndElements();
	}

	size_t threads = environment->OMP_NUM_THREADS;
	eslocal eoffset = _mesh->elements->IDs->datatarray().front();

	std::vector<std::vector<eslocal> > inodes(threads);
	std::vector<std::vector<eslocal> > inodesDistribution(threads);

	inodesDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<eslocal> neighDomains;
		std::vector<eslocal> elements;
		std::vector<eslocal> dnodes;
		std::vector<std::pair<eslocal, eslocal> > intervals;
		const auto &epointers = _mesh->elements->epointers->datatarray();
		for (eslocal d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			neighDomains.clear();
			for (size_t i = 0; i < _mesh->nodes->dintervals[d].size(); i++) {
				auto domains = _mesh->nodes->idomains->cbegin() + _mesh->nodes->dintervals[d][i].pindex;
				neighDomains.insert(neighDomains.end(), domains->begin(), domains->end());
			}
			Esutils::sortAndRemoveDuplicity(neighDomains);

			for (auto nd = neighDomains.begin(); nd != neighDomains.end(); ++nd) {
				if (*nd <= _mesh->elements->firstDomain + d) {
					continue;
				}

				intervals.clear();
				for (size_t i = 0; i < _mesh->nodes->dintervals[d].size(); i++) {
					auto domains = _mesh->nodes->idomains->cbegin() + _mesh->nodes->dintervals[d][i].pindex;
					if (std::binary_search(domains->begin(), domains->end(), *nd)) {
						intervals.push_back(std::make_pair(_mesh->nodes->dintervals[d][i].begin, _mesh->nodes->dintervals[d][i].end));
					}
				}
				for (size_t i = 1; i < intervals.size(); i++) {
					if (intervals[i - 1].second == intervals[i].first) {
						intervals[i].first = intervals[i - 1].first;
						intervals[i - 1].second = intervals[i].second;
					}
				}
				Esutils::removeDuplicity(intervals);

				elements.clear();
				for (size_t i = 0; i < intervals.size(); i++) {
					auto nelements = _mesh->nodes->elements->cbegin() + intervals[i].first;
					for (eslocal e = intervals[i].first; e < intervals[i].second; ++e, ++nelements) {
						for (auto ne = nelements->begin(); ne != nelements->end(); ++ne) {
							if (_mesh->elements->elementsDistribution[d] <= *ne - eoffset && *ne - eoffset < _mesh->elements->elementsDistribution[d + 1]) {
								elements.push_back(*ne - eoffset);
							}
						}
					}
				}
				Esutils::sortAndRemoveDuplicity(elements);

				dnodes.clear();
				for (size_t e = 0; e < elements.size(); ++e) {
					auto nodes = _mesh->elements->nodes->cbegin() + elements[e];
					for (auto f = epointers[elements[e]]->faces->cbegin(); f != epointers[elements[e]]->faces->cend(); ++f) {
						size_t size = 0;
						for (auto fn = f->begin(); fn != f->end(); ++fn) {
							auto it = std::lower_bound(intervals.begin(), intervals.end(), nodes->at(*fn), [&] (const std::pair<eslocal, eslocal> &interval, eslocal node) {
								return interval.second < node;
							});
							if (it != intervals.end() && it->first <= nodes->at(*fn) && nodes->at(*fn) < it->second) {
								++size;
							}
						}
						if (size == f->size()) {
							for (auto fn = f->begin(); fn != f->end(); ++fn) {
								dnodes.push_back(nodes->at(*fn));
							}
						}
					}
				}
				Esutils::sortAndRemoveDuplicity(dnodes);
				inodes[t].insert(inodes[t].end(), dnodes.begin(), dnodes.end());
				inodesDistribution[t].push_back(inodes[t].size());
			}
		}
	}
	Esutils::threadDistributionToFullDistribution(inodesDistribution);

	if (_mesh->sharedInterface != NULL) {
		delete _mesh->sharedInterface;
	}
	_mesh->sharedInterface = new SharedInterfaceStore();
	_mesh->sharedInterface->nodes = new serializededata<eslocal, eslocal>(1, inodes);
	_mesh->sharedInterface->nodeDistribution = inodesDistribution[0];

	finish("computation of shared faces");
}



















