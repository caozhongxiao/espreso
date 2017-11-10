
#include "transformations.h"

#include "../newmesh.h"
#include "../elements/element.h"
#include "../elements/elementstore.h"
#include "../store/domainstore.h"
#include "../store/boundarystore.h"
#include "../store/regionstore.h"

#include "../../basis/point/point.h"
#include "../../basis/containers/serializededata.h"

#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../../wrappers/wmetis.h"
#include "../../wrappers/wparmetis.h"

#include "../../config/ecf/environment.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

void Transformation::reclusterize(NewMesh &mesh)
{
	if (environment->MPIsize == 1) {
		ESINFO(TVERBOSITY) << "Transformation::re-distribution of the mesh to processes skipped (there is only 1 MPI process).";
		return;
	}

	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "Transformation::re-distribution of the mesh to processes started.";

	if (mesh._elems->dual == NULL) {
		Transformation::computeDual(mesh);
	}

	// TODO: ParMetis can get elements coordinates to speed up decomposition
//	if (mesh._elems->coordinates == NULL) {
//		Transformation::computeElementCenters(mesh);
//	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<esglobal> edistribution = mesh._elems->gatherElementDistrubution();
	std::vector<esglobal> partition(mesh._elems->size), permutation(mesh._elems->size), edgeWeights(mesh._elems->dual->datatarray().size());

	size_t edgeConst = 10000;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto dual = mesh._elems->dual->cbegin(t);
		int material;
		Element::TYPE type;

		size_t edgeIndex = mesh._elems->dual->datatarray().distribution()[t];
		for (size_t e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++dual) {
			for (auto neigh = dual->begin(); neigh != dual->end(); ++neigh) {
				auto it = std::lower_bound(mesh._elems->IDs->datatarray().cbegin(), mesh._elems->IDs->datatarray().cend(), *neigh);
				if (it == mesh._elems->IDs->datatarray().cend()) {
					auto halo = std::lower_bound(mesh._halo->IDs->datatarray().cbegin(), mesh._halo->IDs->datatarray().cend(), *neigh);
					material = mesh._halo->material->datatarray()[halo - mesh._halo->IDs->datatarray().cbegin()];
					type = mesh._halo->epointers->datatarray()[halo - mesh._halo->IDs->datatarray().cbegin()]->type;
				} else {
					material = mesh._elems->material->datatarray()[it - mesh._elems->IDs->datatarray().cbegin()];
					type = mesh._elems->epointers->datatarray()[it - mesh._elems->IDs->datatarray().cbegin()]->type;
				}
				edgeWeights[edgeIndex] = 6 * edgeConst + 1;
				if (mesh._elems->epointers->datatarray()[e]->type != type) {
					edgeWeights[edgeIndex] -= 4 * edgeConst;
				}
				if (mesh._elems->material->datatarray()[e] != material) {
					edgeWeights[edgeIndex] -= 2 * edgeConst;
				}
				edgeIndex++;
			}
		}
	}

	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ')<< "Transformation::ParMETIS::KWay started.";
	ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
		edistribution.data(),
		mesh._elems->dual->boundarytaaray().data(), mesh._elems->dual->datatarray().data(),
		0, NULL, // 3, mesh._elems->coordinates->datatarray(),
		0, NULL, edgeWeights.data(),
		partition.data()
	);
	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ')<< "Transformation::ParMETIS::KWay finished.";

	// comment out because weird ParMetis behavior
//	ESINFO(TVERBOSITY) << Info::plain() << "Using ParMETIS to improve edge-cuts: " << edgecut;
//	esglobal prev = 2 * edgecut;
//	while (1.01 * edgecut < prev) {
//		prev = edgecut;
//		edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_AdaptiveRepart,
//			edistribution.data(),
//			mesh._elems->dual->boundarytaaray().data(), mesh._elems->dual->datatarray().data(),
//			0, NULL, // 3, mesh._elems->coordinates->datatarray(),
//			0, NULL, edgeWeights.data(),
//			partition.data()
//		);
//		ESINFO(TVERBOSITY) << Info::plain() << " -> " << edgecut;
//	}
//	ESINFO(TVERBOSITY);

	Transformation::exchangeElements(mesh, partition);
//	Transformation::reindexNodes(mesh);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "Transformation::re-distribution of the mesh to processes finished.";
}

void Transformation::partitiate(NewMesh &mesh, esglobal parts, TFlags::SEPARATE separate)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "Transformation::decomposition of the mesh started.";

	if (mesh._elems->decomposedDual == NULL) {
		Transformation::computeDecomposedDual(mesh, separate);
	}

	size_t threads = environment->OMP_NUM_THREADS;

	auto e2t = [&] (eslocal element) {
		return std::lower_bound(mesh._elems->distribution.begin(), mesh._elems->distribution.end(), element + 1) - mesh._elems->distribution.begin() - 1;
	};

	std::vector<std::vector<int> > part(threads);
	// thread x partID x eID from other part
	std::vector<std::vector<std::vector<eslocal> > > neighElem(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		part[t].resize(mesh._elems->distribution[t + 1] - mesh._elems->distribution[t], -1);
		std::vector<eslocal> stack;
		eslocal current;
		size_t target;
		int partCounter = 0;

		for (size_t i = mesh._elems->distribution[t]; i < mesh._elems->distribution[t + 1]; ++i) {
			if (part[t][i - mesh._elems->distribution[t]] == -1) {
				neighElem[t].push_back(std::vector<eslocal>(threads, -1));
				stack.push_back(i);
				part[t][i - mesh._elems->distribution[t]] = partCounter;
				while (stack.size()) {
					current = stack.back();
					stack.pop_back();
					auto neighs = mesh._elems->decomposedDual->cbegin() + current;
					for (auto e = neighs->begin(); e != neighs->end(); ++e) {
						target = e2t(*e);
						if (target == t) {
							if (part[t][*e - mesh._elems->distribution[t]] == -1) {
								stack.push_back(*e);
								part[t][*e - mesh._elems->distribution[t]] = partCounter;
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
		std::vector<std::pair<size_t, eslocal> > stack; // thread x link to element
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
						int npart = part[current.first][current.second - mesh._elems->distribution[current.first]];
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


	std::vector<eslocal> partition(mesh._elems->size);
	if (nextID == 1) {

		std::vector<eslocal> edgeWeights(mesh._elems->decomposedDual->datatarray().size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto dual = mesh._elems->decomposedDual->cbegin(t);
			int material;
			Element::TYPE type;

			size_t edgeIndex = mesh._elems->decomposedDual->datatarray().distribution()[t];
			for (size_t e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++dual) {
				for (auto neigh = dual->begin(); neigh != dual->end(); ++neigh) {
					auto it = std::lower_bound(mesh._elems->IDs->datatarray().cbegin(), mesh._elems->IDs->datatarray().cend(), *neigh);
					material = mesh._elems->material->datatarray()[it - mesh._elems->IDs->datatarray().cbegin()];
					type = mesh._elems->epointers->datatarray()[it - mesh._elems->IDs->datatarray().cbegin()]->type;
					edgeWeights[edgeIndex] = 6 * edgeConst + 1;

					if (!(separate & TFlags::SEPARATE::ETYPES)) { // dual is already decomposed if separate::ETYPES is true
						if (mesh._elems->epointers->datatarray()[e]->type != type) {
							edgeWeights[edgeIndex] -= 4 * edgeConst;
						}
					}
					if (!(separate & TFlags::SEPARATE::MATERIALS)) { // dual is already decomposed if separate::MATERIALS is true
						if (mesh._elems->material->datatarray()[e] != material) {
							edgeWeights[edgeIndex] -= 2 * edgeConst;
						}
					}
					edgeIndex++;
				}
			}
		}

		ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "Transformation::METIS::KWay started.";
		METIS::call(
				mesh._elems->size,
				mesh._elems->decomposedDual->boundarytaaray().data(), mesh._elems->decomposedDual->datatarray().data(),
				0, NULL, edgeWeights.data(),
				parts, partition.data());
		ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "Transformation::METIS::KWay finished.";
		mesh._domains->clusters.resize(parts, 0);

	} else { // non-continuous dual graph
		// thread x part x elements
		std::vector<std::vector<std::vector<eslocal> > > tdecomposition(threads, std::vector<std::vector<eslocal> >(nextID));
		std::vector<std::vector<eslocal> > tdualsize(threads, std::vector<eslocal>(nextID));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			const auto &dual = mesh._elems->decomposedDual->boundarytaaray();
			for (size_t i = 0, e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++i) {
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
			auto dual = mesh._elems->decomposedDual->cbegin(t);
			size_t partindex;
			int material;
			Element::TYPE type;
			std::vector<eslocal> foffset(nextID), noffset(nextID), edgeIndices(nextID);
			for (int p = 0; p < nextID; p++) {
				foffset[p] = foffsets[p][t];
				edgeIndices[p] = noffset[p] = noffsets[p][t];
			}

			for (size_t i = 0, e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++i, ++dual) {
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

					auto it = std::lower_bound(mesh._elems->IDs->datatarray().cbegin(), mesh._elems->IDs->datatarray().cend(), *node);

					material = mesh._elems->material->datatarray()[it - mesh._elems->IDs->datatarray().cbegin()];
					type = mesh._elems->epointers->datatarray()[it - mesh._elems->IDs->datatarray().cbegin()]->type;
					edgeWeights[partindex][edgeIndices[partindex]] = 6 * edgeConst + 1;

					if (!(separate & TFlags::SEPARATE::ETYPES)) { // dual is already decomposed if separate::ETYPES is true
						if (mesh._elems->epointers->datatarray()[e]->type != type) {
							edgeWeights[partindex][edgeIndices[partindex]] -= 4 * edgeConst;
						}
					}
					if (!(separate & TFlags::SEPARATE::MATERIALS)) { // dual is already decomposed if separate::MATERIALS is true
						if (mesh._elems->material->datatarray()[e] != material) {
							edgeWeights[partindex][edgeIndices[partindex]] -= 2 * edgeConst;
						}
					}

					edgeIndices[partindex]++;
				}
			}
		}

		std::vector<eslocal> pparts(nextID);

		double averageDomainSize = mesh._elems->size / (double)parts;
		size_t partsCounter = 0;
		for (int p = 0; p < nextID; p++) {
			partsCounter += pparts[p] = std::ceil((frames[p].size() - 1) / averageDomainSize);
			mesh._domains->clusters.resize(partsCounter, p);
		}

		ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "Transformation::METIS::KWay started.";
		#pragma omp parallel for
		for (int p = 0; p < nextID; p++) {
			METIS::call(
					frames[p].size() - 1,
					frames[p].data(), neighbors[p].data(),
					0, NULL, edgeWeights[p].data(),
					pparts[p], partition.data() + partoffset[p]);
		}
		ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "Transformation::METIS::KWay finished.";

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

	std::vector<eslocal> domainDistribution;
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
			double averageThreadSize = mesh._elems->size / (double)threads;
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

	Transformation::permuteElements(mesh, permutation, tdistribution);

	std::vector<eslocal> domainCounter(threads);
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

	mesh._domains->size = Esutils::sizesToOffsets(domainCounter);
	domainCounter.push_back(mesh._domains->size);
	mesh._domains->domainDistribution = domainCounter;

	mesh._domains->domainElementBoundaries.push_back(0);
	for (size_t t = 0; t < threads; t++) {
		if (domainDistribution.size() < threads + 1) {
			if (t < domainDistribution.size() - 1) {
				mesh._domains->domainElementBoundaries.push_back(domainDistribution[t + 1]);
			}
		} else {
			auto begin = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t]);
			auto end   = std::lower_bound(domainDistribution.begin(), domainDistribution.end(), tdistribution[t + 1]);
			for (auto it = begin; it != end; ++it) {
				mesh._domains->domainElementBoundaries.push_back(*(it + 1));
			}
		}
	}

	Transformation::arrangeNodes(mesh);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "Transformation::decomposition of the mesh finished.";
}

void Transformation::permuteElements(NewMesh &mesh, const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "Transformation::permutation of elements started.";

	if (mesh._nodes->elems == NULL) {
		Transformation::addLinkFromTo(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
	}

	std::vector<eslocal> backpermutation(permutation.size());
	std::iota(backpermutation.begin(), backpermutation.end(), 0);
	std::sort(backpermutation.begin(), backpermutation.end(), [&] (eslocal i, eslocal j) { return permutation[i] < permutation[j]; });

	size_t threads = environment->OMP_NUM_THREADS;

	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(mesh._neighbours.begin(), mesh._neighbours.end(), neighbor) - mesh._neighbours.begin();
	};

	std::vector<esglobal> IDBoundaries = mesh._elems->gatherElementDistrubution();
	std::vector<std::vector<std::pair<esglobal, esglobal> > > rHalo(mesh._neighbours.size());

	if (mesh._elems->dual != NULL || mesh._nodes->elems != NULL) {
		// thread x neighbor x elements(oldID, newID)
		std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > sHalo(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(mesh._neighbours.size()));

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto ranks = mesh._nodes->ranks->cbegin(t);
			auto elements = mesh._nodes->elems->cbegin(t);
			esglobal begine = IDBoundaries[environment->MPIrank];
			esglobal ende   = IDBoundaries[environment->MPIrank + 1];

			for (auto n = mesh._nodes->distribution[t]; n < mesh._nodes->distribution[t + 1]; ++n, ++ranks, ++elements) {
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

		if (!Communication::exchangeUnknownSize(sHalo[0], rHalo, mesh._neighbours)) {
			ESINFO(ERROR) << "ESPRESO internal error: exchange halo element new IDs while element permutation.";
		}
	}

	auto globalremap = [&] (serializededata<eslocal, esglobal>* data, bool sort) {
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

	auto localremap = [&] (serializededata<eslocal, esglobal>* data, bool sort) {
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

	esglobal firstID = mesh._elems->IDs->datatarray().front();
	mesh._elems->permute(permutation, &distribution);
	std::iota(mesh._elems->IDs->datatarray().begin(), mesh._elems->IDs->datatarray().end(), firstID);

	globalremap(mesh._elems->dual, true);
	globalremap(mesh._nodes->elems, true);
	localremap(mesh._elems->decomposedDual, true);
	localremap(mesh._processBoundaries->elems, false);

	#pragma omp parallel for
	for (size_t i = 0; i < mesh._processBoundaries->elemsIntervals.size(); i++) {
		std::sort(
				mesh._processBoundaries->elems->datatarray().data() + mesh._processBoundaries->elemsIntervals[i].begin,
				mesh._processBoundaries->elems->datatarray().data() + mesh._processBoundaries->elemsIntervals[i].end);
	}

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "Transformation::permutation of elements finished.";
}

void Transformation::exchangeElements(NewMesh &mesh, const std::vector<esglobal> &partition)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "Transformation::exchanging elements started.";

	// 0. Compute targets
	// 1. Serialize element data
	// 2. Serialize node data
	// 3. Send serialized data to target (new MPI processes)
	// 4. Deserialize data
	// 5. Balance nodes data to threads
	// 6. Re-index elements (IDs have to be always increasing)

	if (mesh._nodes->elems == NULL) {
		// need for correctly update nodes ranks
		Transformation::addLinkFromTo(mesh, TFlags::ELEVEL::NODE, TFlags::ELEVEL::ELEMENT);
	}

	size_t threads = environment->OMP_NUM_THREADS;

	// Step 0: Compute targets

	std::vector<int> targets;
	{ // get target processes
		std::vector<std::vector<eslocal> > ttargets(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {

			for (size_t e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e) {
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

	ElementStore *elements = new ElementStore(mesh._eclasses);

	std::vector<std::vector<esglobal> >    elemsIDs(threads);
	std::vector<std::vector<int> >         elemsBody(threads);
	std::vector<std::vector<int> >         elemsMaterial(threads);
	std::vector<std::vector<Element*> > elemsEpointer(threads);
	std::vector<std::vector<eslocal> >     elemsNodesDistribution(threads);
	std::vector<std::vector<esglobal> >    elemsNodesData(threads);
	std::vector<std::vector<int> >         elemsRegions(threads);

	ElementStore *nodes = new ElementStore(mesh._eclasses);

	std::vector<std::vector<esglobal> > nodesIDs(threads);
	std::vector<std::vector<Point> >    nodesCoordinates(threads);
	std::vector<std::vector<eslocal> >  nodesElemsDistribution(threads);
	std::vector<std::vector<esglobal> > nodesElemsData(threads);
	std::vector<std::vector<int> >      nodesRegions(threads);

	// regions are transfered via mask
	int regionsBitMaskSize = mesh._regions.size() / (8 * sizeof(int)) + (mesh._regions.size() % (8 * sizeof(int)) ? 1 : 0);

	// serialize data that have to be exchanged
	// the first thread value denotes the thread data size

	// threads x target x elements(id, body, material, code, dualsize, dualdata, nodesize, nodeindices)
	std::vector<std::vector<std::vector<esglobal> > > sElements(threads, std::vector<std::vector<esglobal> >(targets.size(), std::vector<esglobal>({ 0 })));
	std::vector<std::vector<esglobal> > rElements;

	// threads x target x nodes(id, point, linksize, links)
	std::vector<std::vector<std::vector<esglobal> > > sNodes(threads, std::vector<std::vector<esglobal> >(targets.size(), std::vector<esglobal>({ 0 })));
	std::vector<std::vector<esglobal> > rNodes;

	// Step 1: Serialize element data

	std::vector<int> regionElementMask(mesh._nodes->size * regionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		int maskOffset = 0;
		for (size_t r = 0; r < mesh._regions.size(); r++) {
			maskOffset = r / (8 * sizeof(int));
			if (mesh._regions[r]->etype == TFlags::ELEVEL::ELEMENT) {
				auto begin = std::lower_bound(mesh._regions[r]->elems->datatarray().begin(), mesh._regions[r]->elems->datatarray().end(), mesh._elems->distribution[t]);
				auto end = std::lower_bound(mesh._regions[r]->elems->datatarray().begin(), mesh._regions[r]->elems->datatarray().end(), mesh._elems->distribution[t + 1]);
				for (auto i = begin; i != end; ++i) {
					regionElementMask[*i * regionsBitMaskSize + maskOffset] |= 1 << r;
				}
			}
		}
	}

	elemsNodesDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto IDs = mesh._elems->IDs->datatarray().data();
		auto body = mesh._elems->body->datatarray().data();
		auto material = mesh._elems->material->datatarray().data();
		auto code = mesh._elems->epointers->datatarray().data();
		auto enodes = mesh._elems->nodes->cbegin(t);
		auto nIDs = mesh._nodes->IDs->datatarray().data();

		size_t target;
		for (size_t e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++enodes) {
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
					*(elemsNodesData[t].end() - enodes->size() + n) = mesh._nodes->IDs->datatarray().data()[(*enodes)[n]];
				}
				elemsRegions[t].insert(elemsRegions[t].end(), regionElementMask.begin() + e * regionsBitMaskSize, regionElementMask.begin() + (e + 1) * regionsBitMaskSize);
			} else {
				target = t2i(partition[e]);
				sElements[t][target].insert(sElements[t][target].end(), { IDs[e], body[e], material[e], static_cast<esglobal>(code[e] - mesh._eclasses[t]) });
				sElements[t][target].push_back(enodes->size());
				sElements[t][target].insert(sElements[t][target].end(), enodes->begin(), enodes->end());
				for (size_t n = 0; n < enodes->size(); n++) {
					*(sElements[t][target].end() - enodes->size() + n) = nIDs[(*enodes)[n]];
				}
				sElements[t][target].insert(sElements[t][target].end(), regionElementMask.begin() + e * regionsBitMaskSize, regionElementMask.begin() + (e + 1) * regionsBitMaskSize);
			}
		}
	}

	// Step 2: Serialize node data

	std::vector<int> regionNodeMask(mesh._nodes->size * regionsBitMaskSize);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		int maskOffset = 0;
		for (size_t r = 0; r < mesh._regions.size(); r++) {
			maskOffset = r / (8 * sizeof(int));
			if (mesh._regions[r]->etype == TFlags::ELEVEL::NODE) {
				auto begin = std::lower_bound(mesh._regions[r]->nodes->datatarray().begin(), mesh._regions[r]->nodes->datatarray().end(), mesh._nodes->distribution[t]);
				auto end = std::lower_bound(mesh._regions[r]->nodes->datatarray().begin(), mesh._regions[r]->nodes->datatarray().end(), mesh._nodes->distribution[t + 1]);
				for (auto i = begin; i != end; ++i) {
					regionNodeMask[*i * regionsBitMaskSize + maskOffset] |= 1 << r;
				}
			}
		}
	}

	nodesElemsDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = mesh._nodes->IDs->datatarray();
		const auto &coordinates = mesh._nodes->coordinates->datatarray();
		auto ranks = mesh._nodes->ranks->cbegin(t);
		auto elems = mesh._nodes->elems->cbegin(t);

		const auto &eIDs = mesh._elems->IDs->datatarray();

		size_t target;
		std::vector<bool> last(targets.size() + 1); // targets + me
		for (size_t n = mesh._nodes->distribution[t]; n < mesh._nodes->distribution[t + 1]; ++n, ++elems, ++ranks) {
			std::fill(last.begin(), last.end(), false);
			for (auto e = elems->begin(); e != elems->end(); ++e) {
				auto it = std::lower_bound(eIDs.begin(), eIDs.end(), *e);
				if (it != eIDs.end() && *it == *e) {
					target = t2i(partition[it - eIDs.begin()]);
					if (!last[target] && partition[it - eIDs.begin()] != environment->MPIrank) {
						sNodes[t][target].push_back(IDs[n]);
						sNodes[t][target].insert(sNodes[t][target].end(), sizeof(Point) / sizeof(esglobal), 0);
						memcpy(sNodes[t][target].data() + sNodes[t][target].size() - (sizeof(Point) / sizeof(esglobal)), coordinates.data() + n, sizeof(Point));
						sNodes[t][target].push_back(elems->size());
						sNodes[t][target].insert(sNodes[t][target].end(), elems->begin(), elems->end());
						sNodes[t][target].insert(sNodes[t][target].end(), regionNodeMask.begin() + n * regionsBitMaskSize, regionNodeMask.begin() + (n + 1) * regionsBitMaskSize);
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
						nodesRegions[t].insert(nodesRegions[t].end(), regionNodeMask.begin() + n * regionsBitMaskSize, regionNodeMask.begin() + (n + 1) * regionsBitMaskSize);
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
				elemsEpointer[t].push_back(mesh._eclasses[t] + rElements[i][++e]);
				elemsNodesDistribution[t].push_back(rElements[i][++e]);
				elemsNodesData[t].insert(elemsNodesData[t].end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				e += elemsNodesDistribution[t].back();
				if (elemsNodesDistribution[t].size() > 1) {
					elemsNodesDistribution[t].back() += *(elemsNodesDistribution[t].end() - 2);
				}
				e += regionsBitMaskSize;
				elemsRegions[t].insert(elemsRegions[t].end(), rElements[i].begin() + e * regionsBitMaskSize, rElements[i].begin() + (e + 1) * regionsBitMaskSize);
			}
		}
	}

	// Step 4: Deserialize node data
	std::vector<esglobal> nodeset;
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
		std::vector<std::vector<esglobal> > tnodeset(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			Point point;
			for (size_t n = rdistribution[t]; n + 1 < rdistribution[t + 1]; ) {
				if (!std::binary_search(nodeset.begin(), nodeset.end(), rNodes[i][++n])) {
					nodesIDs[t].push_back(rNodes[i][n]);
					tnodeset[t].push_back(rNodes[i][n]);
					memcpy(&point, rNodes[i].data() + (++n), sizeof(Point));
					nodesCoordinates[t].push_back(point);
					n += sizeof(Point) / sizeof(esglobal);
					nodesElemsDistribution[t].push_back(rNodes[i][n]);
					nodesElemsData[t].insert(nodesElemsData[t].end(), rNodes[i].begin() + n + 1, rNodes[i].begin() + n + 1 + rNodes[i][n]);
					n += nodesElemsDistribution[t].back();
					if (nodesElemsDistribution[t].size() > 1) {
						nodesElemsDistribution[t].back() += *(nodesElemsDistribution[t].end() - 2);
					}
					n += regionsBitMaskSize;
					nodesRegions[t].insert(nodesRegions[t].end(), rNodes[i].begin() + n * regionsBitMaskSize, rNodes[i].begin() + (n + 1) * regionsBitMaskSize);
				} else {
					n += 1 + sizeof(Point) / sizeof(esglobal); // id, Point
					n += rNodes[i][n]; // elems
					n += regionsBitMaskSize; // regions
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

	elements->IDs = new serializededata<eslocal, esglobal>(1, elemsIDs);
	elements->body = new serializededata<eslocal, int>(1, elemsBody);
	elements->material = new serializededata<eslocal, int>(1, elemsMaterial);
	elements->epointers = new serializededata<eslocal, Element*>(1, elemsEpointer);
	elements->nodes = new serializededata<eslocal, esglobal>(elemsNodesDistribution, elemsNodesData); // global IDs

	elements->size = elements->IDs->structures();
	elements->distribution = elements->IDs->datatarray().distribution();

	// Step 5: Balance node data to threads
	std::vector<size_t> nodeDistribution(threads);
	for (size_t t = 1; t < threads; t++) {
		nodeDistribution[t] = nodeDistribution[t - 1] + nodesIDs[t - 1].size();
	}

	serializededata<eslocal, esglobal>::balance(1, nodesIDs);
	serializededata<eslocal, Point>::balance(1, nodesCoordinates);
	serializededata<eslocal, esglobal>::balance(nodesElemsDistribution, nodesElemsData);

	nodes->IDs = new serializededata<eslocal, esglobal>(1, nodesIDs);
	nodes->coordinates = new serializededata<eslocal, Point>(1, nodesCoordinates);
	nodes->elems = new serializededata<eslocal, esglobal>(nodesElemsDistribution, nodesElemsData);
	nodes->size = nodes->IDs->datatarray().size();
	nodes->distribution = nodes->IDs->datatarray().distribution();

	for (size_t r = 0; r < mesh._regions.size(); r++) {
		int maskOffset = r / (8 * sizeof(int));
		if (mesh._regions[r]->etype == TFlags::ELEVEL::NODE) {
			delete mesh._regions[r]->nodes;
			std::vector<std::vector<eslocal> > regionnodes(threads);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = 0; i < nodesRegions[t].size(); i += regionsBitMaskSize) {
					if (nodesRegions[t][i + maskOffset] & (1 << r)) {
						regionnodes[t].push_back(nodeDistribution[t] + i / regionsBitMaskSize);
					}
				}
			}

			serializededata<eslocal, eslocal>::balance(1, regionnodes);
			mesh._regions[r]->nodes = new serializededata<eslocal, eslocal>(1, regionnodes);
			std::sort(mesh._regions[r]->nodes->datatarray().begin(), mesh._regions[r]->nodes->datatarray().end());
		}
		if (mesh._regions[r]->etype == TFlags::ELEVEL::ELEMENT) {
			delete mesh._regions[r]->elems;
			std::vector<std::vector<eslocal> > regionelems(threads);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t i = 0; i < elemsRegions[t].size(); i += regionsBitMaskSize) {
					if (nodesRegions[t][i + maskOffset] & (1 << r)) {
						regionelems[t].push_back(elemDistribution[t] + i / regionsBitMaskSize);
					}
				}
			}

			serializededata<eslocal, eslocal>::balance(1, regionelems);
			mesh._regions[r]->elems = new serializededata<eslocal, eslocal>(1, regionelems);
			std::sort(mesh._regions[r]->elems->datatarray().begin(), mesh._regions[r]->elems->datatarray().end());
		}
	}

	// Step 6: Re-index elements
	// Elements IDs are always kept increasing

	std::vector<esglobal> oldIDBoundaries = mesh._elems->gatherElementDistrubution();
	std::vector<esglobal> newIDBoundaries = elements->gatherElementDistrubution();

	// thread x neighbor x data(ID, new rank)
	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > sHaloTarget(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(mesh._neighbours.size()));
	std::vector<std::vector<std::pair<esglobal, esglobal> > > haloTarget(mesh._neighbours.size());

	// halo elements that will be halo after exchange
	std::vector<std::vector<esglobal> > keptHaloTarget(threads);

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(mesh._neighbours.begin(), mesh._neighbours.end(), neighbour) - mesh._neighbours.begin();
	};

	#pragma omp parallel for
	for (size_t t = 0; t < threads; ++t) {
		auto links = mesh._nodes->elems->cbegin(t);
		auto ranks = mesh._nodes->ranks->cbegin(t);
		std::vector<std::pair<esglobal, esglobal> > buffer;

		for (size_t n = mesh._nodes->distribution[t]; n < mesh._nodes->distribution[t + 1]; ++n, ++links, ++ranks) {
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

	if (!Communication::exchangeUnknownSize(sHaloTarget[0], haloTarget, mesh._neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange elements new MPI processes.";
	}
	Esutils::mergeThreadedUniqueData(haloTarget);

	// haloTarget contains target processes of halo elements
	// fill newIDTargetSBuffet by target processes of elements that was exchanged

	// thread x neighbor x data(oldID, new MPI process)
	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > newIDTargetSBuffet(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(targets.size()));
	std::vector<std::vector<std::pair<esglobal, esglobal> > > newIDTarget(targets.size());

	// elements that are sent to neighbors
	std::vector<std::vector<std::pair<esglobal, esglobal> > > myIDTarget(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto enodes = mesh._elems->nodes->cbegin(t);
		esglobal target;
		bool addToMyIDTargets;

		for (auto e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++enodes) {
			addToMyIDTargets = false;
			if (partition[e] != environment->MPIrank) {
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					auto nelems = mesh._nodes->elems->cbegin() + *n;
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
					myIDTarget[t].push_back(std::make_pair(oldIDBoundaries[environment->MPIrank] + e, partition[e]));
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
	std::sort(newIDTarget[0].begin(), newIDTarget[0].end(), [] (std::pair<esglobal, esglobal> &p1, std::pair<esglobal, esglobal> &p2) {
		if (p1.second == p2.second) {
			return p1.first < p2.first;
		}
		return p1.second < p2.second;
	});

	std::vector<int> neighbors;
	std::vector<std::vector<std::pair<esglobal, esglobal> > > newID, newIDRequest;
	auto begin = newIDTarget[0].begin();
	while (begin != newIDTarget[0].end()) {
		auto end = std::lower_bound(begin, newIDTarget[0].end(), std::make_pair(0, begin->second + 1), [] (std::pair<esglobal, esglobal> &p1, const std::pair<esglobal, esglobal> &p2) { return p1.second < p2.second; });
		if (begin->second != environment->MPIrank) {
			neighbors.push_back(begin->second);
			newID.push_back(std::move(std::vector<std::pair<esglobal, esglobal> >(begin, end)));
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
			auto it = std::lower_bound(permutation.begin(), permutation.end(), newIDRequest[n][i].first, [&] (esglobal i, esglobal val) {
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
		for (auto elem = nodes->elems->begin(t); elem != nodes->elems->end(t); ++elem) {
			for (auto e = elem->begin(); e != elem->end(); ++e) {
				auto it = std::lower_bound(newID[0].begin(), newID[0].end(), std::make_pair(*e, 0));
				if (it == newID[0].end() || it->first != *e) {
					*e = newIDBoundaries[environment->MPIrank] + *std::lower_bound(permutation.begin(), permutation.end(), *e, [&] (esglobal i, esglobal val) {
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
			*n = *std::lower_bound(permutation.begin(), permutation.end(), *n, [&] (esglobal i, esglobal val) {
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
		for (auto elem = nodes->elems->begin(t); elem != nodes->elems->end(t); ++elem) {
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
	std::swap(mesh._elems, elements);
	std::swap(mesh._nodes, nodes);
	delete mesh._halo;
	mesh._halo = new ElementStore(mesh._eclasses);
	mesh._neighbours = neighbors;

	delete elements;
	delete nodes;

	Transformation::computeDual(mesh);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "Transformation::exchanging elements finished.";
}

