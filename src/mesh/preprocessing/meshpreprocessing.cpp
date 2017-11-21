
#include "meshpreprocessing.h"

#include "../mesh.h"

#include "../store/elementstore.h"
#include "../store/nodestore.h"
#include "../store/elementsregionstore.h"
#include "../store/boundaryregionstore.h"

#include "../elements/element.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/logging/logging.h"

#include "../../wrappers/wparmetis.h"
#include "../../wrappers/wmetis.h"

#include <algorithm>
#include <numeric>

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

void MeshPreprocessing::reclusterize()
{
	if (environment->MPIsize == 1) {
		skip("re-distribution of the mesh to processes");
		return;
	}

	start("re-distribution of the mesh to processes");

	if (_mesh->elements->dual == NULL) {
		this->computeDual();
	}

	// TODO: ParMetis can get elements coordinates to speed up decomposition
//	if (_mesh->elements->coordinates == NULL) {
//		Mesh preprocessing :: computeElementCenters(mesh);
//	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<eslocal> edistribution = _mesh->elements->gatherElementsProcDistribution();
	std::vector<eslocal> partition(_mesh->elements->size), permutation(_mesh->elements->size), edgeWeights(_mesh->elements->dual->datatarray().size());

	size_t edgeConst = 10000;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto dual = _mesh->elements->dual->cbegin(t);
		int material;
		Element::TYPE type;

		size_t edgeIndex = _mesh->elements->dual->datatarray().distribution()[t];
		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++dual) {
			for (auto neigh = dual->begin(); neigh != dual->end(); ++neigh) {
				auto it = std::lower_bound(_mesh->elements->IDs->datatarray().cbegin(), _mesh->elements->IDs->datatarray().cend(), *neigh);
				if (it == _mesh->elements->IDs->datatarray().cend()) {
					auto halo = std::lower_bound(_mesh->halo->IDs->datatarray().cbegin(), _mesh->halo->IDs->datatarray().cend(), *neigh);
					material = _mesh->halo->material->datatarray()[halo - _mesh->halo->IDs->datatarray().cbegin()];
					type = _mesh->halo->epointers->datatarray()[halo - _mesh->halo->IDs->datatarray().cbegin()]->type;
				} else {
					material = _mesh->elements->material->datatarray()[it - _mesh->elements->IDs->datatarray().cbegin()];
					type = _mesh->elements->epointers->datatarray()[it - _mesh->elements->IDs->datatarray().cbegin()]->type;
				}
				edgeWeights[edgeIndex] = 6 * edgeConst + 1;
				if (_mesh->elements->epointers->datatarray()[e]->type != type) {
					edgeWeights[edgeIndex] -= 4 * edgeConst;
				}
				if (_mesh->elements->material->datatarray()[e] != material) {
					edgeWeights[edgeIndex] -= 2 * edgeConst;
				}
				edgeIndex++;
			}
		}
	}

	start("ParMETIS::KWay");
	ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
		edistribution.data(),
		_mesh->elements->dual->boundarytaaray().data(), _mesh->elements->dual->datatarray().data(),
		0, NULL, // 3, _mesh->elements->coordinates->datatarray(),
		0, NULL, edgeWeights.data(),
		partition.data()
	);
	finish("ParMETIS::KWay");

	// comment out because weird ParMetis behavior
//	ESINFO(TVERBOSITY) << Info::plain() << "Using ParMETIS to improve edge-cuts: " << edgecut;
//	eslocal prev = 2 * edgecut;
//	while (1.01 * edgecut < prev) {
//		prev = edgecut;
//		edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_AdaptiveRepart,
//			edistribution.data(),
//			_mesh->elements->dual->boundarytaaray().data(), _mesh->elements->dual->datatarray().data(),
//			0, NULL, // 3, _mesh->elements->coordinates->datatarray(),
//			0, NULL, edgeWeights.data(),
//			partition.data()
//		);
//		ESINFO(TVERBOSITY) << Info::plain() << " -> " << edgecut;
//	}
//	ESINFO(TVERBOSITY);

	this->exchangeElements(partition);
//	Mesh preprocessing :: reindexNodes(mesh);

	finish("re-distribution of the mesh to processes");
}

void MeshPreprocessing::partitiate(eslocal parts, bool separateMaterials, bool separateEtype)
{
	start("decomposition of the mesh");

	if (_mesh->elements->decomposedDual == NULL) {
		this->computeDecomposedDual(separateMaterials, separateEtype);
	}

	size_t threads = environment->OMP_NUM_THREADS;

	auto e2t = [&] (eslocal element) ->eslocal {
		return std::lower_bound(_mesh->elements->distribution.begin(), _mesh->elements->distribution.end(), element + 1) - _mesh->elements->distribution.begin() - 1;
	};

	std::vector<std::vector<int> > part(threads);
	// thread x partID x eID from other part
	std::vector<std::vector<std::vector<eslocal> > > neighElem(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		part[t].resize(_mesh->elements->distribution[t + 1] - _mesh->elements->distribution[t], -1);
		std::vector<eslocal> stack;
		eslocal current;
		size_t target;
		int partCounter = 0;

		for (size_t i = _mesh->elements->distribution[t]; i < _mesh->elements->distribution[t + 1]; ++i) {
			if (part[t][i - _mesh->elements->distribution[t]] == -1) {
				neighElem[t].push_back(std::vector<eslocal>(threads, -1));
				stack.push_back(i);
				part[t][i - _mesh->elements->distribution[t]] = partCounter;
				while (stack.size()) {
					current = stack.back();
					stack.pop_back();
					auto neighs = _mesh->elements->decomposedDual->cbegin() + current;
					for (auto e = neighs->begin(); e != neighs->end(); ++e) {
						target = e2t(*e);
						if (target == t) {
							if (part[t][*e - _mesh->elements->distribution[t]] == -1) {
								stack.push_back(*e);
								part[t][*e - _mesh->elements->distribution[t]] = partCounter;
							}
						} else {
							if (neighElem[t][partCounter][target] == -1) {
								neighElem[t][partCounter][target] = *e;
							}
						}
					}
				}
				partCounter++;
			}
		}
	}

	std::vector<std::vector<int> > partID(threads);
	int nextID = 0;
	{ // get parts together
		auto reindexPart = [&] (int oldIndex, int newIndex) {
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = 0; i < partID[t].size(); i++) {
					if (partID[t][i] == oldIndex) {
						partID[t][i] = newIndex;
					}
				}
			}
		};

		for (size_t t = 0; t < threads; t++) {
			partID[t].resize(neighElem[t].size(), -1);
		}
		int reindex;
		std::vector<std::pair<eslocal, eslocal> > stack; // thread x link to element
		std::pair<size_t, eslocal> current;
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = 0; i < neighElem[t].size(); i++) {
				if (partID[t][i] == -1) {
					reindex = -1;
					partID[t][i] = nextID;
					for (size_t n = 0; n < neighElem[t][i].size(); n++) {
						if (neighElem[t][i][n] != -1) {
							stack.push_back(std::make_pair(e2t(neighElem[t][i][n]), neighElem[t][i][n]));
						}
					}
					while (stack.size()) {
						current = stack.back();
						stack.pop_back();
						int npart = part[current.first][current.second - _mesh->elements->distribution[current.first]];
						if (partID[current.first][npart] == -1) {
							partID[current.first][npart] = nextID;
							for (size_t n = 0; n < neighElem[current.first][npart].size(); n++) {
								if (neighElem[current.first][npart][n] != -1) {
									stack.push_back(std::make_pair(e2t(neighElem[current.first][npart][n]), neighElem[current.first][npart][n]));
								}
							}
						}
						if (partID[current.first][npart] < nextID) {
							reindex = partID[current.first][npart];
						}
					}
					if (reindex != -1) {
						reindexPart(nextID, reindex);
					} else {
						++nextID;
					}
				}
			}
		}
	}

	size_t edgeConst = 10000;

	std::vector<int> clusters;
	std::vector<eslocal> partition(_mesh->elements->size);
	if (nextID == 1) {

		std::vector<eslocal> edgeWeights(_mesh->elements->decomposedDual->datatarray().size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto dual = _mesh->elements->decomposedDual->cbegin(t);
			int material;
			Element::TYPE type;

			size_t edgeIndex = _mesh->elements->decomposedDual->datatarray().distribution()[t];
			for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++dual) {
				for (auto neigh = dual->begin(); neigh != dual->end(); ++neigh) {
					auto it = std::lower_bound(_mesh->elements->IDs->datatarray().cbegin(), _mesh->elements->IDs->datatarray().cend(), *neigh);
					material = _mesh->elements->material->datatarray()[it - _mesh->elements->IDs->datatarray().cbegin()];
					type = _mesh->elements->epointers->datatarray()[it - _mesh->elements->IDs->datatarray().cbegin()]->type;
					edgeWeights[edgeIndex] = 6 * edgeConst + 1;

					if (!separateEtype) { // dual is already decomposed if separate::ETYPES is true
						if (_mesh->elements->epointers->datatarray()[e]->type != type) {
							edgeWeights[edgeIndex] -= 4 * edgeConst;
						}
					}
					if (!separateMaterials) { // dual is already decomposed if separate::MATERIALS is true
						if (_mesh->elements->material->datatarray()[e] != material) {
							edgeWeights[edgeIndex] -= 2 * edgeConst;
						}
					}
					edgeIndex++;
				}
			}
		}

		start("METIS::KWay");
		METIS::call(
				_mesh->elements->size,
				_mesh->elements->decomposedDual->boundarytaaray().data(), _mesh->elements->decomposedDual->datatarray().data(),
				0, NULL, edgeWeights.data(),
				parts, partition.data());
		finish("METIS::KWay");
		clusters.resize(parts, 0);

	} else { // non-continuous dual graph
		// thread x part x elements
		std::vector<std::vector<std::vector<eslocal> > > tdecomposition(threads, std::vector<std::vector<eslocal> >(nextID));
		std::vector<std::vector<eslocal> > tdualsize(threads, std::vector<eslocal>(nextID));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			const auto &dual = _mesh->elements->decomposedDual->boundarytaaray();
			for (size_t i = 0, e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++i) {
				tdecomposition[t][partID[t][part[t][i]]].push_back(e);
				tdualsize[t][partID[t][part[t][i]]] += dual[e + 1] - dual[e];
			}
		}
		std::vector<std::vector<eslocal> > foffsets(nextID), noffsets(nextID);
		std::vector<eslocal> partoffset(nextID);
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			foffsets[p].push_back(0);
			noffsets[p].push_back(0);
			for (size_t t = 1; t < threads; t++) {
				foffsets[p].push_back(tdecomposition[0][p].size());
				noffsets[p].push_back(tdualsize[0][p]);
				tdecomposition[0][p].insert(tdecomposition[0][p].end(), tdecomposition[t][p].begin(), tdecomposition[t][p].end());
				tdualsize[0][p] += tdualsize[t][p];
			}
		}
		for (int p = 1; p < nextID; p++) {
			partoffset[p] = partoffset[p - 1] + tdecomposition[0][p - 1].size();
		}

		std::vector<std::vector<eslocal> > frames(nextID), neighbors(nextID), edgeWeights(nextID);
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			frames[p].resize(1 + tdecomposition[0][p].size());
			neighbors[p].resize(tdualsize[0][p]);
			edgeWeights[p].resize(tdualsize[0][p]);
		}

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto dual = _mesh->elements->decomposedDual->cbegin(t);
			size_t partindex;
			int material;
			Element::TYPE type;
			std::vector<eslocal> foffset(nextID), noffset(nextID), edgeIndices(nextID);
			for (int p = 0; p < nextID; p++) {
				foffset[p] = foffsets[p][t];
				edgeIndices[p] = noffset[p] = noffsets[p][t];
			}

			for (size_t i = 0, e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++i, ++dual) {
				partindex = partID[t][part[t][i]];

				frames[partindex][++foffset[partindex]] = dual->size();
				if (i) {
					frames[partindex][foffset[partindex]] += frames[partindex][foffset[partindex] - 1];
				} else {
					frames[partindex][foffset[partindex]] += noffset[partindex];
				}
				auto node = dual->begin();
				for (eslocal n = frames[partindex][foffset[partindex]] - dual->size(); n < frames[partindex][foffset[partindex]]; ++n, ++node) {
					neighbors[partindex][n] = std::lower_bound(tdecomposition[0][partindex].begin(), tdecomposition[0][partindex].end(), *node) - tdecomposition[0][partindex].begin();

					auto it = std::lower_bound(_mesh->elements->IDs->datatarray().cbegin(), _mesh->elements->IDs->datatarray().cend(), *node);

					material = _mesh->elements->material->datatarray()[it - _mesh->elements->IDs->datatarray().cbegin()];
					type = _mesh->elements->epointers->datatarray()[it - _mesh->elements->IDs->datatarray().cbegin()]->type;
					edgeWeights[partindex][edgeIndices[partindex]] = 6 * edgeConst + 1;

					if (!separateEtype) { // dual is already decomposed if separate::ETYPES is true
						if (_mesh->elements->epointers->datatarray()[e]->type != type) {
							edgeWeights[partindex][edgeIndices[partindex]] -= 4 * edgeConst;
						}
					}
					if (!separateMaterials) { // dual is already decomposed if separate::MATERIALS is true
						if (_mesh->elements->material->datatarray()[e] != material) {
							edgeWeights[partindex][edgeIndices[partindex]] -= 2 * edgeConst;
						}
					}

					edgeIndices[partindex]++;
				}
			}
		}

		std::vector<eslocal> pparts(nextID);

		double averageDomainSize = _mesh->elements->size / (double)parts;
		size_t partsCounter = 0;
		for (int p = 0; p < nextID; p++) {
			partsCounter += pparts[p] = std::ceil((frames[p].size() - 1) / averageDomainSize);
			clusters.resize(partsCounter, p);
		}

		start("METIS::KWay");
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			METIS::call(
					frames[p].size() - 1,
					frames[p].data(), neighbors[p].data(),
					0, NULL, edgeWeights[p].data(),
					pparts[p], partition.data() + partoffset[p]);
		}
		finish("METIS::KWay");

		std::vector<eslocal> ppartition = partition;
		nextID = 0;
		for (size_t p = 0; p < tdecomposition[0].size(); p++) {
			for (size_t i = 0; i < tdecomposition[0][p].size(); ++i) {
				partition[tdecomposition[0][p][i]] = ppartition[partoffset[p] + i] + nextID;
			}
			nextID += pparts[p];
		}
	}

	std::vector<eslocal> permutation(partition.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		if (partition[i] == partition[j]) {
			return i < j;
		}
		return partition[i] < partition[j];
	});

	std::vector<size_t> domainDistribution;
	std::vector<size_t> tdistribution;

	eslocal partindex = 0;
	auto begin = permutation.begin();
	while (begin != permutation.end()) {
		domainDistribution.push_back(begin - permutation.begin());
		begin = std::lower_bound(begin, permutation.end(), ++partindex, [&] (eslocal i, eslocal val) {
			return partition[i] < val;
		});
	}
	domainDistribution.push_back(permutation.size());

	// TODO: improve domain distribution for more complicated decomposition
	if (domainDistribution.size() == threads + 1) {
		tdistribution = std::vector<size_t>(domainDistribution.begin(), domainDistribution.end());
	} else {
		if (domainDistribution.size() < threads + 1) {
			tdistribution = tarray<eslocal>::distribute(threads, permutation.size());
		} else {
			double averageThreadSize = _mesh->elements->size / (double)threads;
			tdistribution.push_back(0);
			for (size_t t = 0; t < threads - 1; t++) {
				auto more = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution.back() + averageThreadSize);
				auto less = more - 1;
				if (std::fabs(*less - averageThreadSize * (t + 1)) < std::fabs(*more - averageThreadSize * (t + 1))) {
					tdistribution.push_back(*less);
				} else {
					tdistribution.push_back(*more);
				}
			}
			tdistribution.push_back(permutation.size());
		}
	}

	this->permuteElements(permutation, tdistribution);

	for (size_t i = 1; i < domainDistribution.size(); i++) {
		if (domainDistribution[i - 1] != domainDistribution[i]) {
			_mesh->elements->clusters.push_back(clusters[i - 1]);
		}
	}
	Esutils::removeDuplicity(domainDistribution);

	std::vector<size_t> domainCounter(threads);
	for (size_t t = 0; t < threads; t++) {
		if (domainDistribution.size() < threads + 1) {
			if (t < domainDistribution.size() - 1) {
				++domainCounter[t];
			}
		} else {
			auto begin = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t + 1]);
			for (auto it = begin; it != end; ++it) {
				++domainCounter[t];
			}
		}
	}

	_mesh->elements->ndomains = Esutils::sizesToOffsets(domainCounter);
	_mesh->elements->firstDomain = _mesh->elements->ndomains;
	Communication::exscan(_mesh->elements->firstDomain);
	domainCounter.push_back(_mesh->elements->ndomains);
	_mesh->elements->domainDistribution = domainCounter;

	_mesh->elements->elementsDistribution.push_back(0);
	for (size_t t = 0; t < threads; t++) {
		if (domainDistribution.size() < threads + 1) {
			if (t < domainDistribution.size() - 1) {
				_mesh->elements->elementsDistribution.push_back(domainDistribution[t + 1]);
			}
		} else {
			auto begin = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t + 1]);
			for (auto it = begin; it != end; ++it) {
				_mesh->elements->elementsDistribution.push_back(*(it + 1));
			}
		}
	}

	MeshPreprocessing::arrangeNodes();

	finish("decomposition of the mesh");
}

void MeshPreprocessing::exchangeElements(const std::vector<eslocal> &partition)
{
	start("exchanging elements");

	// 0. Compute targets
	// 1. Serialize element data
	// 2. Serialize node data
	// 3. Send serialized data to target (new MPI processes)
	// 4. Deserialize data
	// 5. Balance nodes data to threads
	// 6. Re-index elements (IDs have to be always increasing)

	if (_mesh->nodes->elements == NULL) {
		// need for correctly update nodes ranks
		this->linkNodesAndElements();
	}

	size_t threads = environment->OMP_NUM_THREADS;

	// Step 0: Compute targets

	std::vector<int> targets;
	{ // get target processes
		std::vector<std::vector<eslocal> > ttargets(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e) {
				if (partition[e] != environment->MPIrank) {
					// assumes small number of targets
					auto it = std::lower_bound(ttargets[t].begin(), ttargets[t].end(), partition[e]);
					if (it == ttargets[t].end() || *it != partition[e]) {
						ttargets[t].insert(it, partition[e]);
					}
				}
			}
		}

		Esutils::mergeThreadedUniqueData(ttargets);
		targets = ttargets[0];
	}

	auto t2i = [ & ] (size_t target) {
		return std::lower_bound(targets.begin(), targets.end(), target) - targets.begin();
	};

	ElementStore *elements = new ElementStore(_mesh->_eclasses);

	std::vector<std::vector<eslocal> >    elemsIDs(threads);
	std::vector<std::vector<int> >         elemsBody(threads);
	std::vector<std::vector<int> >         elemsMaterial(threads);
	std::vector<std::vector<Element*> > elemsEpointer(threads);
	std::vector<std::vector<eslocal> >     elemsNodesDistribution(threads);
	std::vector<std::vector<eslocal> >    elemsNodesData(threads);
	std::vector<std::vector<int> >         elemsRegions(threads);

	NodeStore *nodes = new NodeStore();

	std::vector<std::vector<eslocal> >  nodesIDs(threads);
	std::vector<std::vector<Point> >    nodesCoordinates(threads);
	std::vector<std::vector<eslocal> >  nodesElemsDistribution(threads);
	std::vector<std::vector<eslocal> >  nodesElemsData(threads);
	std::vector<std::vector<int> >      nodesRegions(threads);

	// regions are transfered via mask
	int eregionsBitMaskSize = _mesh->elementsRegions.size() / (8 * sizeof(int)) + (_mesh->elementsRegions.size() % (8 * sizeof(int)) ? 1 : 0);
	int bregionsBitMaskSize = _mesh->boundaryRegions.size() / (8 * sizeof(int)) + (_mesh->boundaryRegions.size() % (8 * sizeof(int)) ? 1 : 0);

	// serialize data that have to be exchanged
	// the first thread value denotes the thread data size

	// threads x target x elements(id, body, material, code, dualsize, dualdata, nodesize, nodeindices)
	std::vector<std::vector<std::vector<eslocal> > > sElements(threads, std::vector<std::vector<eslocal> >(targets.size(), std::vector<eslocal>({ 0 })));
	std::vector<std::vector<eslocal> > rElements;

	// threads x target x nodes(id, point, linksize, links)
	std::vector<std::vector<std::vector<eslocal> > > sNodes(threads, std::vector<std::vector<eslocal> >(targets.size(), std::vector<eslocal>({ 0 })));
	std::vector<std::vector<eslocal> > rNodes;

	// Step 1: Serialize element data

	std::vector<int> regionElementMask(_mesh->nodes->size * eregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		int maskOffset = 0;
		for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(int));
			auto begin = std::lower_bound(_mesh->elementsRegions[r]->elements->datatarray().begin(), _mesh->elementsRegions[r]->elements->datatarray().end(), _mesh->elements->distribution[t]);
			auto end = std::lower_bound(_mesh->elementsRegions[r]->elements->datatarray().begin(), _mesh->elementsRegions[r]->elements->datatarray().end(), _mesh->elements->distribution[t + 1]);
			for (auto i = begin; i != end; ++i) {
				regionElementMask[*i * eregionsBitMaskSize + maskOffset] |= 1 << r;
			}
		}
	}

	elemsNodesDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto IDs = _mesh->elements->IDs->datatarray().data();
		auto body = _mesh->elements->body->datatarray().data();
		auto material = _mesh->elements->material->datatarray().data();
		auto code = _mesh->elements->epointers->datatarray().data();
		auto enodes = _mesh->elements->nodes->cbegin(t);
		auto nIDs = _mesh->nodes->IDs->datatarray().data();

		size_t target;
		for (size_t e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++enodes) {
			if (partition[e] == environment->MPIrank) {
				elemsIDs[t].push_back(IDs[e]);
				elemsBody[t].push_back(body[e]);
				elemsMaterial[t].push_back(material[e]);
				elemsEpointer[t].push_back(code[e]);
				elemsNodesDistribution[t].push_back(enodes->size());
				elemsNodesData[t].insert(elemsNodesData[t].end(), enodes->begin(), enodes->end());
				if (elemsNodesDistribution[t].size() > 1) {
					elemsNodesDistribution[t].back() += *(elemsNodesDistribution[t].end() - 2);
				}
				for (size_t n = 0; n < enodes->size(); n++) {
					*(elemsNodesData[t].end() - enodes->size() + n) = _mesh->nodes->IDs->datatarray().data()[(*enodes)[n]];
				}
				elemsRegions[t].insert(elemsRegions[t].end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			} else {
				target = t2i(partition[e]);
				sElements[t][target].insert(sElements[t][target].end(), { IDs[e], body[e], material[e], static_cast<eslocal>(code[e] - _mesh->_eclasses[t]) });
				sElements[t][target].push_back(enodes->size());
				sElements[t][target].insert(sElements[t][target].end(), enodes->begin(), enodes->end());
				for (size_t n = 0; n < enodes->size(); n++) {
					*(sElements[t][target].end() - enodes->size() + n) = nIDs[(*enodes)[n]];
				}
				sElements[t][target].insert(sElements[t][target].end(), regionElementMask.begin() + e * eregionsBitMaskSize, regionElementMask.begin() + (e + 1) * eregionsBitMaskSize);
			}
		}
	}

	// Step 2: Serialize node data

	std::vector<int> regionNodeMask(_mesh->nodes->size * bregionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		int maskOffset = 0;
		for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
			maskOffset = r / (8 * sizeof(int));
			auto begin = std::lower_bound(_mesh->boundaryRegions[r]->nodes->datatarray().begin(), _mesh->boundaryRegions[r]->nodes->datatarray().end(), _mesh->nodes->distribution[t]);
			auto end = std::lower_bound(_mesh->boundaryRegions[r]->nodes->datatarray().begin(), _mesh->boundaryRegions[r]->nodes->datatarray().end(), _mesh->nodes->distribution[t + 1]);
			for (auto i = begin; i != end; ++i) {
				regionNodeMask[*i * bregionsBitMaskSize + maskOffset] |= 1 << r;
			}
		}
	}

	nodesElemsDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = _mesh->nodes->IDs->datatarray();
		const auto &coordinates = _mesh->nodes->coordinates->datatarray();
		auto ranks = _mesh->nodes->ranks->cbegin(t);
		auto elems = _mesh->nodes->elements->cbegin(t);

		const auto &eIDs = _mesh->elements->IDs->datatarray();

		size_t target;
		std::vector<bool> last(targets.size() + 1); // targets + me
		for (size_t n = _mesh->nodes->distribution[t]; n < _mesh->nodes->distribution[t + 1]; ++n, ++elems, ++ranks) {
			std::fill(last.begin(), last.end(), false);
			for (auto e = elems->begin(); e != elems->end(); ++e) {
				auto it = std::lower_bound(eIDs.begin(), eIDs.end(), *e);
				if (it != eIDs.end() && *it == *e) {
					target = t2i(partition[it - eIDs.begin()]);
					if (!last[target] && partition[it - eIDs.begin()] != environment->MPIrank) {
						sNodes[t][target].push_back(IDs[n]);
						sNodes[t][target].insert(sNodes[t][target].end(), sizeof(Point) / sizeof(eslocal), 0);
						memcpy(sNodes[t][target].data() + sNodes[t][target].size() - (sizeof(Point) / sizeof(eslocal)), coordinates.data() + n, sizeof(Point));
						sNodes[t][target].push_back(elems->size());
						sNodes[t][target].insert(sNodes[t][target].end(), elems->begin(), elems->end());
						sNodes[t][target].insert(sNodes[t][target].end(), regionNodeMask.begin() + n * bregionsBitMaskSize, regionNodeMask.begin() + (n + 1) * bregionsBitMaskSize);
						last[target] = true;
					}
					if (!last.back() && partition[it - eIDs.begin()] == environment->MPIrank) {
						nodesIDs[t].push_back(IDs[n]);
						nodesCoordinates[t].push_back(coordinates[n]);
						nodesElemsDistribution[t].push_back(elems->size());
						nodesElemsData[t].insert(nodesElemsData[t].end(), elems->begin(), elems->end());
						if (nodesElemsDistribution[t].size() > 1) {
							nodesElemsDistribution[t].back() += *(nodesElemsDistribution[t].end() - 2);
						}
						nodesRegions[t].insert(nodesRegions[t].end(), regionNodeMask.begin() + n * bregionsBitMaskSize, regionNodeMask.begin() + (n + 1) * bregionsBitMaskSize);
						last.back() = true;
					}
				}
			}
		}
	}

	// Step 3: Send data to target processes

	#pragma omp parallel for
	for (size_t target = 0; target < targets.size(); ++target) {
		sElements[0][target].front() = sElements[0][target].size();
		sNodes[0][target].front() = sNodes[0][target].size();
		for (size_t t = 1; t < threads; t++) {
			sElements[t][target].front() = sElements[t][target].size();
			sElements[0][target].insert(sElements[0][target].end(), sElements[t][target].begin(), sElements[t][target].end());
			sNodes[t][target].front() = sNodes[t][target].size();
			sNodes[0][target].insert(sNodes[0][target].end(), sNodes[t][target].begin(), sNodes[t][target].end());
		}
	}

	if (!Communication::sendVariousTargets(sElements[0], rElements, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange elements data.";
	}

	if (!Communication::sendVariousTargets(sNodes[0], rNodes, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange nodes data.";
	}

	// Step 4: Deserialize element data
	for (size_t i = 0; i < rElements.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rElements[i].size()) {
			rdistribution.push_back(p += rElements[i][p]);
		}
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = rdistribution[t]; e + 1 < rdistribution[t + 1]; ) {
				elemsIDs[t].push_back(rElements[i][++e]);
				elemsBody[t].push_back(rElements[i][++e]);
				elemsMaterial[t].push_back(rElements[i][++e]);
				elemsEpointer[t].push_back(_mesh->_eclasses[t] + rElements[i][++e]);
				elemsNodesDistribution[t].push_back(rElements[i][++e]);
				elemsNodesData[t].insert(elemsNodesData[t].end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				e += elemsNodesDistribution[t].back();
				if (elemsNodesDistribution[t].size() > 1) {
					elemsNodesDistribution[t].back() += *(elemsNodesDistribution[t].end() - 2);
				}
				e += eregionsBitMaskSize;
				elemsRegions[t].insert(elemsRegions[t].end(), rElements[i].begin() + e * eregionsBitMaskSize, rElements[i].begin() + (e + 1) * eregionsBitMaskSize);
			}
		}
	}

	// Step 4: Deserialize node data
	std::vector<eslocal> nodeset;
	for (size_t t = 0; t < threads; t++) {
		nodeset.insert(nodeset.end(), nodesIDs[t].begin(), nodesIDs[t].end());
	}
	std::sort(nodeset.begin(), nodeset.end());
	for (size_t i = 0; i < rNodes.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rNodes[i].size()) {
			rdistribution.push_back(p += rNodes[i][p]);
		}
		std::vector<std::vector<eslocal> > tnodeset(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			Point point;
			for (size_t n = rdistribution[t]; n + 1 < rdistribution[t + 1]; ) {
				if (!std::binary_search(nodeset.begin(), nodeset.end(), rNodes[i][++n])) {
					nodesIDs[t].push_back(rNodes[i][n]);
					tnodeset[t].push_back(rNodes[i][n]);
					memcpy(&point, rNodes[i].data() + (++n), sizeof(Point));
					nodesCoordinates[t].push_back(point);
					n += sizeof(Point) / sizeof(eslocal);
					nodesElemsDistribution[t].push_back(rNodes[i][n]);
					nodesElemsData[t].insert(nodesElemsData[t].end(), rNodes[i].begin() + n + 1, rNodes[i].begin() + n + 1 + rNodes[i][n]);
					n += nodesElemsDistribution[t].back();
					if (nodesElemsDistribution[t].size() > 1) {
						nodesElemsDistribution[t].back() += *(nodesElemsDistribution[t].end() - 2);
					}
					n += bregionsBitMaskSize;
					nodesRegions[t].insert(nodesRegions[t].end(), rNodes[i].begin() + n * bregionsBitMaskSize, rNodes[i].begin() + (n + 1) * bregionsBitMaskSize);
				} else {
					n += 1 + sizeof(Point) / sizeof(eslocal); // id, Point
					n += rNodes[i][n]; // elems
					n += bregionsBitMaskSize; // regions
				}
			}
		}
		for (size_t t = 0; t < threads; t++) {
			nodeset.insert(nodeset.end(), tnodeset[t].begin(), tnodeset[t].end());
		}
		std::sort(nodeset.begin(), nodeset.end());
	}

	Esutils::threadDistributionToFullDistribution(elemsNodesDistribution);
	Esutils::threadDistributionToFullDistribution(nodesElemsDistribution);

	// elements are redistributed later while decomposition -> distribution is not changed now
	std::vector<size_t> elemDistribution(threads);
	for (size_t t = 1; t < threads; t++) {
		elemDistribution[t] = elemDistribution[t - 1] + elemsIDs[t - 1].size();
	}

	elements->IDs = new serializededata<eslocal, eslocal>(1, elemsIDs);
	elements->body = new serializededata<eslocal, int>(1, elemsBody);
	elements->material = new serializededata<eslocal, int>(1, elemsMaterial);
	elements->epointers = new serializededata<eslocal, Element*>(1, elemsEpointer);
	elements->nodes = new serializededata<eslocal, eslocal>(elemsNodesDistribution, elemsNodesData); // global IDs

	elements->size = elements->IDs->structures();
	elements->distribution = elements->IDs->datatarray().distribution();

	// Step 5: Balance node data to threads
	std::vector<size_t> nodeDistribution(threads);
	for (size_t t = 1; t < threads; t++) {
		nodeDistribution[t] = nodeDistribution[t - 1] + nodesIDs[t - 1].size();
	}

	serializededata<eslocal, eslocal>::balance(1, nodesIDs);
	serializededata<eslocal, Point>::balance(1, nodesCoordinates);
	serializededata<eslocal, eslocal>::balance(nodesElemsDistribution, nodesElemsData);

	nodes->IDs = new serializededata<eslocal, eslocal>(1, nodesIDs);
	nodes->coordinates = new serializededata<eslocal, Point>(1, nodesCoordinates);
	nodes->elements = new serializededata<eslocal, eslocal>(nodesElemsDistribution, nodesElemsData);
	nodes->size = nodes->IDs->datatarray().size();
	nodes->distribution = nodes->IDs->datatarray().distribution();

	for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
		int maskOffset = r / (8 * sizeof(int));
		delete _mesh->elementsRegions[r]->elements;
		std::vector<std::vector<eslocal> > regionelems(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = 0; i < elemsRegions[t].size(); i += eregionsBitMaskSize) {
				if (nodesRegions[t][i + maskOffset] & (1 << r)) {
					regionelems[t].push_back(elemDistribution[t] + i / eregionsBitMaskSize);
				}
			}
		}

		serializededata<eslocal, eslocal>::balance(1, regionelems);
		_mesh->elementsRegions[r]->elements = new serializededata<eslocal, eslocal>(1, regionelems);
		std::sort(_mesh->elementsRegions[r]->elements->datatarray().begin(), _mesh->elementsRegions[r]->elements->datatarray().end());
	}

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		int maskOffset = r / (8 * sizeof(int));
		delete _mesh->boundaryRegions[r]->nodes;
		std::vector<std::vector<eslocal> > regionnodes(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = 0; i < nodesRegions[t].size(); i += bregionsBitMaskSize) {
				if (nodesRegions[t][i + maskOffset] & (1 << r)) {
					regionnodes[t].push_back(nodeDistribution[t] + i / bregionsBitMaskSize);
				}
			}
		}

		serializededata<eslocal, eslocal>::balance(1, regionnodes);
		_mesh->boundaryRegions[r]->nodes = new serializededata<eslocal, eslocal>(1, regionnodes);
		std::sort(_mesh->boundaryRegions[r]->nodes->datatarray().begin(), _mesh->boundaryRegions[r]->nodes->datatarray().end());
	}

	// Step 6: Re-index elements
	// Elements IDs are always kept increasing

	std::vector<eslocal> oldIDBoundaries = _mesh->elements->gatherElementsProcDistribution();
	std::vector<eslocal> newIDBoundaries = elements->gatherElementsProcDistribution();

	// thread x neighbor x data(ID, new rank)
	std::vector<std::vector<std::vector<std::pair<eslocal, eslocal> > > > sHaloTarget(threads, std::vector<std::vector<std::pair<eslocal, eslocal> > >(_mesh->neighbours.size()));
	std::vector<std::vector<std::pair<eslocal, eslocal> > > haloTarget(_mesh->neighbours.size());

	// halo elements that will be halo after exchange
	std::vector<std::vector<eslocal> > keptHaloTarget(threads);

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbour) - _mesh->neighbours.begin();
	};

	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		auto links = _mesh->nodes->elements->cbegin(t);
		auto ranks = _mesh->nodes->ranks->cbegin(t);
		std::vector<std::pair<eslocal, eslocal> > buffer;

		for (size_t n = _mesh->nodes->distribution[t]; n < _mesh->nodes->distribution[t + 1]; ++n, ++links, ++ranks) {
			if (ranks->size() > 1) {
				buffer.clear();
				bool keep = false;
				size_t ksize = keptHaloTarget[t].size();
				for (auto e = links->begin(); e != links->end(); ++e) {
					if (oldIDBoundaries[environment->MPIrank] <= *e && *e < oldIDBoundaries[environment->MPIrank + 1]) {
						buffer.push_back(std::make_pair(*e, partition[*e - oldIDBoundaries[environment->MPIrank]]));
						if (partition[*e - oldIDBoundaries[environment->MPIrank]] == environment->MPIrank) {
							keep = true;
						}
					} else {
						keptHaloTarget[t].push_back(*e);
					}
				}
				if (!keep) {
					keptHaloTarget[t].resize(ksize);
				}
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != environment->MPIrank) {
						sHaloTarget[t][n2i(*rank)].insert(sHaloTarget[t][n2i(*rank)].end(), buffer.begin(), buffer.end());
					}
				}
			}
		}

		Esutils::sortAndRemoveDuplicity(sHaloTarget[t]);
		Esutils::sortAndRemoveDuplicity(keptHaloTarget[t]);
	}
	Esutils::mergeThreadedUniqueData(sHaloTarget);

	if (!Communication::exchangeUnknownSize(sHaloTarget[0], haloTarget, _mesh->neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange elements new MPI processes.";
	}
	Esutils::mergeThreadedUniqueData(haloTarget);

	// haloTarget contains target processes of halo elements
	// fill newIDTargetSBuffet by target processes of elements that was exchanged

	// thread x neighbor x data(oldID, new MPI process)
	std::vector<std::vector<std::vector<std::pair<eslocal, eslocal> > > > newIDTargetSBuffet(threads, std::vector<std::vector<std::pair<eslocal, eslocal> > >(targets.size()));
	std::vector<std::vector<std::pair<eslocal, eslocal> > > newIDTarget(targets.size());

	// elements that are sent to neighbors
	std::vector<std::vector<std::pair<eslocal, eslocal> > > myIDTarget(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto enodes = _mesh->elements->nodes->cbegin(t);
		eslocal target;
		bool addToMyIDTargets;

		for (auto e = _mesh->elements->distribution[t]; e < _mesh->elements->distribution[t + 1]; ++e, ++enodes) {
			addToMyIDTargets = false;
			if (partition[e] != environment->MPIrank) {
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					auto nelems = _mesh->nodes->elements->cbegin() + *n;
					for (auto ne = nelems->begin(); ne != nelems->end(); ++ne) {
						if (oldIDBoundaries[environment->MPIrank] <= *ne && *ne < oldIDBoundaries[environment->MPIrank + 1]) {
							target = partition[*ne - oldIDBoundaries[environment->MPIrank]];
						} else {
							target = std::lower_bound(haloTarget[0].begin(), haloTarget[0].end(), std::make_pair(*ne, 0))->second;
						}
						if (target != partition[e]) {
							newIDTargetSBuffet[t][t2i(partition[e])].push_back(std::make_pair(*ne, target));
						}
						if (target == environment->MPIrank) {
							addToMyIDTargets = true;
						}
					}
				}
				if (addToMyIDTargets) {
					myIDTarget[t].push_back(std::make_pair((eslocal)(oldIDBoundaries[environment->MPIrank] + e), partition[e]));
				}
			}
		}
		Esutils::sortAndRemoveDuplicity(newIDTargetSBuffet[t]);
		Esutils::sortAndRemoveDuplicity(myIDTarget[t]);
	}

	Esutils::mergeThreadedUniqueData(newIDTargetSBuffet);

	if (!Communication::sendVariousTargets(newIDTargetSBuffet[0], newIDTarget, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange elements new IDs.";
	}
	if (newIDTarget.size() < threads) {
		newIDTarget.resize(threads);
	}
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = 0; i < keptHaloTarget[t].size(); i++) {
			auto it = std::lower_bound(haloTarget[0].begin(), haloTarget[0].end(), std::make_pair(keptHaloTarget[t][i], 0));
			newIDTarget[t].push_back(*it);
		}
	}
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		newIDTarget[t].insert(newIDTarget[t].end(), myIDTarget[t].begin(), myIDTarget[t].end());
	}

	Esutils::mergeThreadedUniqueData(newIDTarget);

	// Ask targets to new IDs
	std::sort(newIDTarget[0].begin(), newIDTarget[0].end(), [] (std::pair<eslocal, eslocal> &p1, std::pair<eslocal, eslocal> &p2) {
		if (p1.second == p2.second) {
			return p1.first < p2.first;
		}
		return p1.second < p2.second;
	});

	std::vector<int> neighbors;
	std::vector<std::vector<std::pair<eslocal, eslocal> > > newID, newIDRequest;
	auto begin = newIDTarget[0].begin();
	while (begin != newIDTarget[0].end()) {
		auto end = std::lower_bound(begin, newIDTarget[0].end(), std::make_pair(0, begin->second + 1), [] (std::pair<eslocal, eslocal> &p1, const std::pair<eslocal, eslocal> &p2) { return p1.second < p2.second; });
		if (begin->second != environment->MPIrank) {
			neighbors.push_back(begin->second);
			newID.push_back(std::move(std::vector<std::pair<eslocal, eslocal> >(begin, end)));
		}
		begin = end;
	}

	newIDRequest.resize(newID.size());
	if (!Communication::exchangeUnknownSize(newID, newIDRequest, neighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: get new ID requests";
	}

	std::vector<eslocal> permutation(elements->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return elements->IDs->datatarray().data()[i] < elements->IDs->datatarray().data()[j]; });

	#pragma omp parallel for
	for (size_t n = 0; n < neighbors.size(); n++) {
		for (size_t i = 0; i < newIDRequest[n].size(); i++) {
			auto it = std::lower_bound(permutation.begin(), permutation.end(), newIDRequest[n][i].first, [&] (eslocal i, eslocal val) {
				return elements->IDs->datatarray().data()[i] < val;
			});
			if (it == permutation.end()) {
				ESINFO(ERROR) << "ESPRESO internal error: request for unknown element ID.";
			}
			newIDRequest[n][i].second = newIDBoundaries[environment->MPIrank] + *it;
		}
	}

	if (!Communication::exchangeKnownSize(newIDRequest, newID, neighbors)) {
		ESINFO(ERROR) << "ESPRESO internal error: get new elements IDs.";
	}
	Esutils::mergeThreadedUniqueData(newID);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto elem = nodes->elements->begin(t); elem != nodes->elements->end(t); ++elem) {
			for (auto e = elem->begin(); e != elem->end(); ++e) {
				auto it = std::lower_bound(newID[0].begin(), newID[0].end(), std::make_pair(*e, 0));
				if (it == newID[0].end() || it->first != *e) {
					*e = newIDBoundaries[environment->MPIrank] + *std::lower_bound(permutation.begin(), permutation.end(), *e, [&] (eslocal i, eslocal val) {
						return elements->IDs->datatarray().data()[i] < val;
					});
				} else {
					*e = it->second;
				}
			}
			std::sort(elem->begin(), elem->end());
		}
	}

	permutation.resize(nodes->size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return nodes->IDs->datatarray().data()[i] < nodes->IDs->datatarray().data()[j]; });

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = elements->nodes->begin(t)->begin(); n != elements->nodes->end(t)->begin(); ++n) {
			*n = *std::lower_bound(permutation.begin(), permutation.end(), *n, [&] (eslocal i, eslocal val) {
				return nodes->IDs->datatarray().data()[i] < val;
			});
		}
	}

	std::vector<std::vector<eslocal> > rankBoundaries(threads);
	std::vector<std::vector<int> > rankData(threads);

	rankBoundaries.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		int rank, rsize;
		for (auto elem = nodes->elements->begin(t); elem != nodes->elements->end(t); ++elem) {
			rsize = 1;
			rankData[t].push_back(std::lower_bound(newIDBoundaries.begin(), newIDBoundaries.end(), *elem->begin() + 1) - newIDBoundaries.begin() - 1);
			for (auto e = elem->begin() + 1; e != elem->end(); ++e) {
				rank = std::lower_bound(newIDBoundaries.begin(), newIDBoundaries.end(), *e + 1) - newIDBoundaries.begin() - 1;
				if (rank != rankData[t].back()) {
					rankData[t].push_back(rank);
					++rsize;
				}
			}
			rankBoundaries[t].push_back(rsize);
			if (rankBoundaries[t].size() > 1) {
				rankBoundaries[t].back() += *(rankBoundaries[t].end() - 2);
			}
		}
	}

	Esutils::threadDistributionToFullDistribution(rankBoundaries);

	nodes->ranks = new serializededata<eslocal, int>(rankBoundaries, rankData);

	std::iota(elements->IDs->datatarray().begin(), elements->IDs->datatarray().end(), newIDBoundaries[environment->MPIrank]);
	std::swap(_mesh->elements, elements);
	std::swap(_mesh->nodes, nodes);
	delete _mesh->halo;
	_mesh->halo = new ElementStore(_mesh->_eclasses);
	_mesh->neighbours = neighbors;

	delete elements;
	delete nodes;

	this->computeDual();

	finish("exchanging elements");
}

void MeshPreprocessing::permuteElements(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
{
	start("permutation of elements");

	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	std::vector<eslocal> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (eslocal i, eslocal j) { return permutation[i] < permutation[j]; });

	size_t threads = environment->OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(_mesh->neighbours.begin(), _mesh->neighbours.end(), neighbor) - _mesh->neighbours.begin();
	};

	std::vector<eslocal> IDBoundaries = _mesh->elements->gatherElementsProcDistribution();
	std::vector<std::vector<std::pair<eslocal, eslocal> > > rHalo(_mesh->neighbours.size());

	if (_mesh->elements->dual != NULL || _mesh->nodes->elements != NULL) {
		// thread x neighbor x elements(oldID, newID)
		std::vector<std::vector<std::vector<std::pair<eslocal, eslocal> > > > sHalo(threads, std::vector<std::vector<std::pair<eslocal, eslocal> > >(_mesh->neighbours.size()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto ranks = _mesh->nodes->ranks->cbegin(t);
			auto elements = _mesh->nodes->elements->cbegin(t);
			eslocal begine = IDBoundaries[environment->MPIrank];
			eslocal ende   = IDBoundaries[environment->MPIrank + 1];

			for (auto n = _mesh->nodes->distribution[t]; n < _mesh->nodes->distribution[t + 1]; ++n, ++ranks, ++elements) {
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != environment->MPIrank) {
						for (auto e = elements->begin(); e != elements->end(); ++e) {
							if (begine <= *e && *e < ende) {
								sHalo[t][n2i(*rank)].push_back(std::make_pair(*e, backpermutation[*e - begine] + begine));
							}
						}
					}
				}
			}

			for (size_t n = 0; n < sHalo[t].size(); ++n) {
				Esutils::sortAndRemoveDuplicity(sHalo[t][n]);
			}
		}

		Esutils::mergeThreadedUniqueData(sHalo);

		if (!Communication::exchangeUnknownSize(sHalo[0], rHalo, _mesh->neighbours)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange halo element new IDs while element permutation.";
		}
	}

	auto globalremap = [&] (serializededata<eslocal, eslocal>* data, bool sort) {
		if (data == NULL) {
			return;
		}
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			int source;
			for (auto e = data->begin(t); e != data->end(t); ++e) {
				for (auto n = e->begin(); n != e->end(); ++n) {
					source = std::lower_bound(IDBoundaries.begin(), IDBoundaries.end(), *n + 1) - IDBoundaries.begin() - 1;
					if (source == environment->MPIrank) {
						*n = IDBoundaries[environment->MPIrank] + backpermutation[*n - IDBoundaries[environment->MPIrank]];
					} else {
						*n = std::lower_bound(rHalo[n2i(source)].begin(), rHalo[n2i(source)].end(), std::make_pair(*n, 0))->second;
					}
				}
				if (sort) {
					std::sort(e->begin(), e->end());
				}
			}
		}
	};

	auto localremap = [&] (serializededata<eslocal, eslocal>* data, bool sort) {
		if (data == NULL) {
			return;
		}
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (auto e = data->begin(t); e != data->end(t); ++e) {
				for (auto n = e->begin(); n != e->end(); ++n) {
					*n = backpermutation[*n];
				}
				if (sort) {
					std::sort(e->begin(), e->end());
				}
			}
		}
	};

	eslocal firstID = _mesh->elements->IDs->datatarray().front();
	_mesh->elements->permute(permutation, distribution);
	std::iota(_mesh->elements->IDs->datatarray().begin(), _mesh->elements->IDs->datatarray().end(), firstID);

	globalremap(_mesh->elements->dual, true);
	globalremap(_mesh->nodes->elements, true);
	localremap(_mesh->elements->decomposedDual, true);

	// TODO: MESH
//	localremap(_mesh->_processBoundaries->elems, false);
//
//	#pragma omp parallel for
//	for (size_t i = 0; i < _mesh->_processBoundaries->elemsIntervals.size(); i++) {
//		std::sort(
//				_mesh->_processBoundaries->elems->datatarray().data() + _mesh->_processBoundaries->elemsIntervals[i].begin,
//				_mesh->_processBoundaries->elems->datatarray().data() + _mesh->_processBoundaries->elemsIntervals[i].end);
//	}

	finish("permutation of elements");
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

	std::vector<size_t> IDBoundaries = _mesh->elements->gatherElementsDistribution();

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

				auto facepointer = epointer->front()->facepointers->cbegin(t);
				for (auto face = epointer->front()->faces->cbegin(t); face != epointer->front()->faces->cend(t); ++face, ++facepointer) {

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

void MeshPreprocessing::arrangeNodes()
{
	if (_mesh->elements->domainDistribution.size() == 0) {
		skip("arrange nodes");
		return;
	}

	start("arrange nodes");

	if (_mesh->nodes->elements == NULL) {
		this->linkNodesAndElements();
	}

	std::vector<eslocal> externalBoundary, internalBoundary;
	this->computeBoundaryNodes(externalBoundary, internalBoundary);

	size_t threads = environment->OMP_NUM_THREADS;

	// TODO: element nodes reindex to domains
//	std::vector<std::vector<eslocal> > domainNodes(_mesh->elements->ndomains);
//
//	mesh._domains->elems = new serializededata<eslocal, eslocal>(*mesh._elems->nodes);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		for (eslocal d = mesh._domains->domainDistribution[t]; d < mesh._domains->domainDistribution[t + 1]; ++d) {
//			domainNodes[d].insert(domainNodes[d].end(),
//					(mesh._domains->elems->cbegin() + mesh._domains->domainElementBoundaries[d])->begin(),
//					(mesh._domains->elems->cbegin() + mesh._domains->domainElementBoundaries[d + 1])->begin());
//
//			Esutils::sortAndRemoveDuplicity(domainNodes[d]);
//			for (
//				auto n = (mesh._domains->elems->begin() + mesh._domains->domainElementBoundaries[d])->begin();
//				n != (mesh._domains->elems->begin() + mesh._domains->domainElementBoundaries[d + 1])->begin();
//				++n) {
//
//				*n = std::lower_bound(domainNodes[d].begin(), domainNodes[d].end(), *n) - domainNodes[d].begin();
//			}
//		}
//	}
//
//	mesh._domains->domainNodeBoundaries.clear();
//	mesh._domains->domainNodeBoundaries.push_back(0);
//	for (eslocal d = 0; d < mesh._domains->size; ++d) {
//		mesh._domains->domainNodeBoundaries.push_back(mesh._domains->domainNodeBoundaries.back() + domainNodes[d].size());
//	}
//
//	std::vector<std::vector<eslocal> > tdomainNodes(threads);
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		for (eslocal d = mesh._domains->domainDistribution[t]; d < mesh._domains->domainDistribution[t + 1]; ++d) {
//			tdomainNodes[t].insert(tdomainNodes[t].end(), domainNodes[d].begin(), domainNodes[d].end());
//		}
//	}
//
//	mesh._domains->nodes = new serializededata<eslocal, eslocal>(1, tdomainNodes);

	std::vector<size_t> eDist = _mesh->elements->gatherElementsDistribution();
	std::vector<eslocal> dProcDist = _mesh->elements->gatherDomainsProcDistribution();
	std::vector<std::vector<eslocal> > domainsDistribution(threads), domainsData(threads);
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
					domainsProcs[t].push_back(std::lower_bound(dProcDist.begin(), dProcDist.end(), domainsData[t].back() + 1) - dProcDist.begin()  -1);
				}
				prev = domain;
			}

			domainsDistribution[t].push_back(domainsData[t].size());
		}
	}
	Esutils::threadDistributionToFullDistribution(domainsDistribution);

	for (size_t t = 1; t < threads; t++) {
		domainsDistribution[0].insert(domainsDistribution[0].end(), domainsDistribution[t].begin(), domainsDistribution[t].end());
		domainsData[0].insert(domainsData[0].end(), domainsData[t].begin(), domainsData[t].end());
		domainsProcs[0].insert(domainsProcs[0].end(), domainsProcs[t].begin(), domainsProcs[t].end());
	}

	std::vector<eslocal> permutation;
	permutation.reserve(_mesh->nodes->size);
	permutation.insert(permutation.end(), externalBoundary.begin(), externalBoundary.end());
	permutation.insert(permutation.end(), internalBoundary.begin(), internalBoundary.end());

	auto ebpointer = externalBoundary.begin();
	auto ibpointer = internalBoundary.begin();
	for (eslocal n = 0; n < _mesh->nodes->size; ++n) {
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

	auto comp = [&] (eslocal i, eslocal j) {
		eslocal di = domainsDistribution[0][i], isize = domainsDistribution[0][i + 1] - domainsDistribution[0][i];
		eslocal dj = domainsDistribution[0][j], jsize = domainsDistribution[0][j + 1] - domainsDistribution[0][j];

		if (isize == jsize) {
			for (eslocal d = 0; d < isize; d++) {
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

	auto equalNeighs = [&] (eslocal i, eslocal j) {
		eslocal di = domainsDistribution[0][i], isize = domainsDistribution[0][i + 1] - domainsDistribution[0][i];
		eslocal dj = domainsDistribution[0][j], jsize = domainsDistribution[0][j + 1] - domainsDistribution[0][j];
		if (isize != jsize) {
			return false;
		}
		for (size_t n = 0; n < isize; ++n) {
			if (domainsData[0][di + n] != domainsData[0][dj + n]) {
				return false;
			}
		}
		return true;
	};

	std::vector<std::vector<eslocal> > iBoundary(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = _mesh->nodes->distribution[t]; i < _mesh->nodes->distribution[t + 1]; ++i) {
			if (i > 0 && !equalNeighs(permutation[i], permutation[i - 1])) {
				iBoundary[t].push_back(i);
			}
		}
	}
	iBoundary[0].push_back(externalBoundary.size());
	iBoundary[0].push_back(externalBoundary.size() + internalBoundary.size());

	Esutils::mergeThreadedUniqueData(iBoundary);
	if (iBoundary[0].back() == _mesh->nodes->size) {
		iBoundary[0].pop_back();
	}

	std::vector<eslocal> intervalDomainsDistribution, intervalDomainsData;
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
	_mesh->nodes->idomains = new serializededata<eslocal, eslocal>(tarray<eslocal>(0, threads, intervalDomainsDistribution), tarray<eslocal>(0, threads, intervalDomainsData));
	_mesh->nodes->iprocesses = new serializededata<eslocal, int>(tarray<eslocal>(0, threads, intervalDomainsDistribution), tarray<eslocal>(0, threads, intervalDomainsProcs));

	_mesh->nodes->externalIntervals = 0;
	auto iti = _mesh->nodes->pintervals.begin();
	while (iti != _mesh->nodes->pintervals.end() && iti->end <= (eslocal)externalBoundary.size()) {
		++_mesh->nodes->externalIntervals;
		++iti;
	}

	#pragma omp parallel for
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i) {
		std::sort(permutation.begin() + _mesh->nodes->pintervals[i].begin, permutation.begin() + _mesh->nodes->pintervals[i].end, [&] (eslocal i, eslocal j) {
			return _mesh->nodes->IDs->datatarray()[i] < _mesh->nodes->IDs->datatarray()[j];
		});
	}

	// 0 -> external + inner, 1 -> inner, 2 -> external, 3 -> rest
	auto getBoundary = [&] (eslocal i) {
		auto neighbors = _mesh->nodes->idomains->cbegin() + i;
		if (i < _mesh->nodes->externalIntervals) {
			if (neighbors->size() > 1) {
				return 0;
			} else {
				return 2;
			}
		} else {
			if (neighbors->size() > 1) {
				return 1;
			} else {
				return 3;
			}
		}
	};

	std::vector<eslocal> ipermutation(_mesh->nodes->pintervals.size());
	std::iota(ipermutation.begin(), ipermutation.end(), 0);
	std::sort(ipermutation.begin(), ipermutation.end(), [&] (eslocal i, eslocal j) {
		int b1 = getBoundary(i), b2 = getBoundary(j);
		auto n1 = _mesh->nodes->idomains->cbegin() + i;
		auto n2 = _mesh->nodes->idomains->cbegin() + j;
		if (b1 == b2) {
			if (n1->size() != n2->size()) {
				return n1->size() > n2->size();
			}

			for (eslocal d = 0; d < n1->size(); d++) {
				if ((*n1)[d] != (*n2)[d]) {
					return (*n1)[d] < (*n2)[d];
				}
			}
		}
		return b1 < b2;
	});

	_mesh->nodes->idomains->permute(ipermutation, _mesh->nodes->idomains->boundarytaaray().distribution());
	_mesh->nodes->iprocesses->permute(ipermutation, _mesh->nodes->iprocesses->boundarytaaray().distribution());
	std::vector<ProcessInterval> permutedIntervals;
	for (size_t i = 0; i < _mesh->nodes->pintervals.size(); i++) {
		permutedIntervals.push_back(_mesh->nodes->pintervals[ipermutation[i]]);
	}

	std::vector<eslocal> finalpermutation;
	finalpermutation.reserve(permutation.size());

	_mesh->nodes->pintervals.clear();
	ipermutation.clear();
	eslocal ioffset = 0;
	for (size_t t = 0; t < threads; t++) {
		for (size_t d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			auto domains = _mesh->nodes->idomains->cbegin();
			for (size_t i = 0; i < permutedIntervals.size(); ++i, ++domains) {
				if (*std::lower_bound(domains->begin(), domains->end(), _mesh->elements->firstDomain) == _mesh->elements->firstDomain + d) {
					finalpermutation.insert(finalpermutation.end(), permutation.begin() + permutedIntervals[i].begin, permutation.begin() + permutedIntervals[i].end);
					eslocal isize = permutedIntervals[i].end - permutedIntervals[i].begin;
					permutedIntervals[i].begin = ioffset;
					ioffset += isize;
					permutedIntervals[i].end = ioffset;
					ipermutation.push_back(i);
					_mesh->nodes->pintervals.push_back(permutedIntervals[i]);
				}
			}
		}
	}
	_mesh->nodes->idomains->permute(ipermutation, _mesh->nodes->idomains->boundarytaaray().distribution());
	_mesh->nodes->iprocesses->permute(ipermutation, _mesh->nodes->iprocesses->boundarytaaray().distribution());

	_mesh->nodes->dintervals.resize(_mesh->elements->ndomains);
	_mesh->nodes->gintervals.resize(_mesh->elements->ndomains);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
			eslocal doffset = 0;
			auto neighbors = _mesh->nodes->idomains->cbegin();
			for (size_t i = 0; i < _mesh->nodes->pintervals.size(); ++i, ++neighbors) {
				auto iit = std::lower_bound(neighbors->begin(), neighbors->end(), _mesh->elements->firstDomain + d);
				if (iit != neighbors->end() && *iit == _mesh->elements->firstDomain + d) {
					_mesh->nodes->dintervals[d].push_back(DomainInterval(_mesh->nodes->pintervals[i].begin, _mesh->nodes->pintervals[i].end, i, doffset));
					_mesh->nodes->gintervals[d].push_back(GluingInterval(_mesh->nodes->dintervals[d].back(), iit - neighbors->begin(), neighbors->size()));
					doffset += _mesh->nodes->pintervals[i].end - _mesh->nodes->pintervals[i].begin;
				}
			}
		}
	}

	_mesh->nodes->permute(finalpermutation);

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

	localremap(_mesh->elements->nodes);

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		localremap(_mesh->boundaryRegions[r]->nodes);
		std::sort(_mesh->boundaryRegions[r]->nodes->datatarray().begin(), _mesh->boundaryRegions[r]->nodes->datatarray().end());

		{
			auto offset = _mesh->boundaryRegions[r]->nodes->datatarray().cbegin();
			auto cend = _mesh->boundaryRegions[r]->nodes->datatarray().cend();
			for (size_t i = 0; i < _mesh->nodes->pintervals.size(); i++) {
				auto begin = offset = std::lower_bound(offset, cend, _mesh->nodes->pintervals[i].begin);
				auto end   = offset = std::lower_bound(offset, cend, _mesh->nodes->pintervals[i].end);
				_mesh->boundaryRegions[r]->nodesIntervals.push_back(_mesh->nodes->pintervals[i]);
				_mesh->boundaryRegions[r]->nodesIntervals.back().begin = begin - _mesh->boundaryRegions[r]->nodes->datatarray().cbegin();
				_mesh->boundaryRegions[r]->nodesIntervals.back().end   = end   - _mesh->boundaryRegions[r]->nodes->datatarray().cbegin();
			}
		}

		_mesh->boundaryRegions[r]->domainNodesIntervals.resize(_mesh->elements->ndomains);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t d = _mesh->elements->domainDistribution[t]; d < _mesh->elements->domainDistribution[t + 1]; d++) {
				auto offset = _mesh->boundaryRegions[r]->nodes->datatarray().cbegin();
				auto cend = _mesh->boundaryRegions[r]->nodes->datatarray().cend();
				for (size_t i = 0; i < _mesh->nodes->dintervals[d].size(); i++) {
					auto begin = offset = std::lower_bound(offset, cend, _mesh->nodes->dintervals[d][i].begin);
					auto end   = offset = std::lower_bound(offset, cend, _mesh->nodes->dintervals[d][i].end);
					_mesh->boundaryRegions[r]->domainNodesIntervals[d].push_back(_mesh->nodes->dintervals[d][i]);
					_mesh->boundaryRegions[r]->domainNodesIntervals[d].back().begin = begin - _mesh->boundaryRegions[r]->nodes->datatarray().cbegin();
					_mesh->boundaryRegions[r]->domainNodesIntervals[d].back().end   = end   - _mesh->boundaryRegions[r]->nodes->datatarray().cbegin();
				}

			}
		}

	}

	finish("arrange nodes");
}
