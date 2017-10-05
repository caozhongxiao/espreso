
#include "transformations.h"

#include "../newmesh.h"
#include "../elements/elementstore.h"

#include "../../basis/point/point.h"
#include "../../basis/containers/serializededata.h"

#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../../wrappers/wmetis.h"
#include "../../wrappers/wparmetis.h"

#include "../../config/ecf/environment.h"

#include <iostream>
#include <algorithm>
#include <numeric>

using namespace espreso;

static void printElement(ElementStore *store)
{
	for (int rank = 0; rank < environment->MPIsize; rank++) {
		if (rank == environment->MPIrank) {
			std::cout << "RANK: " << rank << "\n";
			for (size_t t = 0; t < store->distribution.size() - 1; t++) {
				// std::cout << "T: " << t << "\n";
				auto nodes = store->nodes->cbegin(t);
				auto IDs = store->IDs->datatarray().data();
				auto dual = store->dual->cbegin(t);
				for (size_t e = store->distribution[t]; e < store->distribution[t + 1]; ++e, ++nodes, ++dual) {
					std::cout << "ID(" << IDs[e] << "): ";
					for (auto n = dual->begin(); n != dual->end(); ++n) {
						std::cout << "ID(" << *n << ") ";
					}
//					for (size_t n = 0; n < nodes->size(); n++) {
//						std::cout << (*nodes)[n] << " ";
//					}
					std::cout << "\n";
				}
			}
		}
		MPI_Barrier(environment->MPICommunicator);
	}
	MPI_Barrier(environment->MPICommunicator);
}

static void printNodes(ElementStore *store)
{
	for (int rank = 0; rank < environment->MPIsize; rank++) {
		if (rank == environment->MPIrank) {
			std::cout << "RANK: " << rank << "\n";
			for (size_t t = 0; t < store->distribution.size() - 1; t++) {
				std::cout << "T: " << t << "\n";
				auto coordinates = store->coordinates->datatarray().data();
				auto IDs = store->IDs->datatarray().data();

				for (size_t e = store->distribution[t]; e < store->distribution[t + 1]; ++e) {
					std::cout << IDs[e] << ": " << coordinates[e] << "\n";
				}
			}
		}
		MPI_Barrier(environment->MPICommunicator);
	}
	MPI_Barrier(environment->MPICommunicator);
}

template <typename TType>
static void printVector(const std::vector<TType> &vector)
{
	for (int rank = 0; rank < environment->MPIsize; rank++) {
		if (rank == environment->MPIrank) {
			std::cout << "RANK: " << rank << "\n";
			std::cout << vector;
		}
		MPI_Barrier(environment->MPICommunicator);
	}
	MPI_Barrier(environment->MPICommunicator);
}


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

	std::vector<esglobal> edistribution(environment->MPIsize + 1);

	esglobal esize = mesh._elems->size;
	Communication::exscan(esize);

	MPI_Allgather(&esize, sizeof(esglobal), MPI_BYTE, edistribution.data(), sizeof(esglobal), MPI_BYTE, MPI_COMM_WORLD);
	edistribution.back() = esize + mesh._elems->size;
	MPI_Bcast(&edistribution.back(), sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);


	std::vector<esglobal> partition(mesh._elems->size), permutation(mesh._elems->size), edgeWeights(mesh._elems->dual->datatarray().size());

	size_t edgeConst = 10000;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto dual = mesh._elems->dual->cbegin(t);
		int material;
		NewElement::TYPE type;

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

	esglobal edgecut = ParMETIS::call(ParMETIS::METHOD::ParMETIS_V3_PartKway,
		edistribution.data(),
		mesh._elems->dual->boundarytaaray().data(), mesh._elems->dual->datatarray().data(),
		0, NULL, // 3, mesh._elems->coordinates->datatarray(),
		0, NULL, edgeWeights.data(),
		partition.data()
	);

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

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "Transformation::re-distribution of the mesh to processes finished.";
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

	ElementStore *elements = new ElementStore();

	std::vector<std::vector<esglobal> >    elemsIDs(threads);
	std::vector<std::vector<int> >         elemsBody(threads);
	std::vector<std::vector<int> >         elemsMaterial(threads);
	std::vector<std::vector<NewElement*> > elemsEpointer(threads);
	std::vector<std::vector<eslocal> >     elemsNodesDistribution(threads);
	std::vector<std::vector<esglobal> >    elemsNodesData(threads);

	ElementStore *nodes = new ElementStore();

	std::vector<std::vector<esglobal> > nodesIDs(threads);
	std::vector<std::vector<Point> >    nodesCoordinates(threads);
	std::vector<std::vector<eslocal> >  nodesElemsDistribution(threads);
	std::vector<std::vector<esglobal> > nodesElemsData(threads);

	// serialize data that have to be exchanged
	// the first thread value denotes the thread data size

	// threads x target x elements(id, body, material, code, dualsize, dualdata, nodesize, nodeindices)
	std::vector<std::vector<std::vector<esglobal> > > sElements(threads, std::vector<std::vector<esglobal> >(targets.size(), std::vector<esglobal>({ 0 })));
	std::vector<std::vector<esglobal> > rElements;

	// threads x target x nodes(id, point, linksize, links)
	std::vector<std::vector<std::vector<esglobal> > > sNodes(threads, std::vector<std::vector<esglobal> >(targets.size(), std::vector<esglobal>({ 0 })));
	std::vector<std::vector<esglobal> > rNodes;

	// Step 1: Serialize element data

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
			} else {
				target = t2i(partition[e]);

				sElements[t][target].insert(sElements[t][target].end(), { IDs[e], body[e], material[e], static_cast<esglobal>(code[e] - mesh._eclasses[t].data()) });
				sElements[t][target].push_back(enodes->size());
				sElements[t][target].insert(sElements[t][target].end(), enodes->begin(), enodes->end());
				for (size_t n = 0; n < enodes->size(); n++) {
					*(sElements[t][target].end() - enodes->size() + n) = nIDs[(*enodes)[n]];
				}
			}
		}
	}

	// Step 2: Serialize node data

	nodesElemsDistribution.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto IDs = mesh._nodes->IDs->datatarray().data();
		auto coordinates = mesh._nodes->coordinates->datatarray().data();
		auto ranks = mesh._nodes->ranks->cbegin(t);
		auto elems = mesh._nodes->elems->cbegin(t);

		auto eIDs = mesh._elems->IDs->datatarray();

		size_t target;
		std::vector<bool> last(targets.size() + 1); // targets + me
		for (size_t n = mesh._nodes->distribution[t]; n < mesh._nodes->distribution[t + 1]; ++n, ++elems, ++ranks) {
			std::fill(last.begin(), last.end(), false);
			for (auto e = elems->begin(); e != elems->end(); ++e) {
				auto it = std::lower_bound(eIDs.begin(), eIDs.end(), *e);
				if (it != eIDs.end() && *it == *e) {
					target = t2i(partition[it - eIDs.begin()]);
					if (!last[target] && !std::binary_search(ranks->begin(), ranks->end(), partition[it - eIDs.begin()])) {
						sNodes[t][target].push_back(IDs[n]);
						sNodes[t][target].insert(sNodes[t][target].end(), sizeof(Point) / sizeof(esglobal), 0);
						memcpy(sNodes[t][target].data() + sNodes[t][target].size() - (sizeof(Point) / sizeof(esglobal)), coordinates + n, sizeof(Point));
						sNodes[t][target].push_back(elems->size());
						sNodes[t][target].insert(sNodes[t][target].end(), elems->begin(), elems->end());
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
				elemsEpointer[t].push_back(mesh._eclasses[t].data() + rElements[i][++e]);
				elemsNodesDistribution[t].push_back(rElements[i][++e]);
				elemsNodesData[t].insert(elemsNodesData[t].end(), rElements[i].begin() + e + 1, rElements[i].begin() + e + 1 + rElements[i][e]);
				e += elemsNodesDistribution[t].back();
				if (elemsNodesDistribution[t].size() > 1) {
					elemsNodesDistribution[t].back() += *(elemsNodesDistribution[t].end() - 2);
				}
			}
		}
	}

	// Step 4: Deserialize node data
	std::vector<esglobal> nodeset;
	for (size_t i = 0; i < rNodes.size(); i++) {
		std::vector<size_t> rdistribution({ 0 });
		size_t p = 0;
		while (p < rNodes[i].size()) {
			rdistribution.push_back(p += rNodes[i][p]);
		}
		std::vector<std::vector<esglobal> > tnodeset(threads);
		Point point;
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
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
				} else {
					n += 1 + sizeof(Point) / sizeof(esglobal); // id, Point
					n += rNodes[i][++n]; // elems
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
	elements->IDs = new serializededata<eslocal, esglobal>(1, elemsIDs);
	elements->body = new serializededata<eslocal, int>(1, elemsBody);
	elements->material = new serializededata<eslocal, int>(1, elemsMaterial);
	elements->nodes = new serializededata<eslocal, esglobal>(elemsNodesDistribution, elemsNodesData); // global IDs

	elements->size = elements->IDs->structures();
	elements->distribution = elements->IDs->datatarray().distribution();


	// Step 5: Balance node data to threads
	serializededata<eslocal, esglobal>::balance(1, nodesIDs);
	serializededata<eslocal, Point>::balance(1, nodesCoordinates);
	serializededata<eslocal, esglobal>::balance(nodesElemsDistribution, nodesElemsData);

	nodes->IDs = new serializededata<eslocal, esglobal>(1, nodesIDs);
	nodes->coordinates = new serializededata<eslocal, Point>(1, nodesCoordinates);
	nodes->elems = new serializededata<eslocal, esglobal>(nodesElemsDistribution, nodesElemsData);
	nodes->size = nodes->IDs->datatarray().size();
	nodes->distribution = nodes->IDs->datatarray().distribution();

	// Step 6: Re-index elements
	// Elements IDs are always kept increasing

	std::vector<esglobal> oldIDBoundaries(environment->MPIsize + 1);
	esglobal esize = mesh._elems->size;
	Communication::exscan(esize);

	MPI_Allgather(&esize, sizeof(esglobal), MPI_BYTE, oldIDBoundaries.data(), sizeof(esglobal), MPI_BYTE, MPI_COMM_WORLD);
	oldIDBoundaries.back() = esize + mesh._elems->size;
	MPI_Bcast(&oldIDBoundaries.back(), sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);

	std::vector<esglobal> newIDBoundaries(environment->MPIsize + 1);
	esize = elements->size;
	Communication::exscan(esize);

	MPI_Allgather(&esize, sizeof(esglobal), MPI_BYTE, newIDBoundaries.data(), sizeof(esglobal), MPI_BYTE, MPI_COMM_WORLD);
	newIDBoundaries.back() = esize + elements->size;
	MPI_Bcast(&newIDBoundaries.back(), sizeof(esglobal), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);

	// thread x neighbor x data(ID, new rank)
	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > sHaloTarget(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(mesh._neighbours.size()));
	std::vector<std::vector<std::pair<esglobal, esglobal> > > haloTarget(mesh._neighbours.size());

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
				for (auto e = links->begin(); e != links->end(); ++e) {
					if (oldIDBoundaries[environment->MPIrank] <= *e && *e < oldIDBoundaries[environment->MPIrank + 1]) {
						buffer.push_back(std::make_pair(*e, partition[*e - oldIDBoundaries[environment->MPIrank]]));
					}
				}
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != environment->MPIrank) {
						sHaloTarget[t][n2i(*rank)].insert(sHaloTarget[t][n2i(*rank)].end(), buffer.begin(), buffer.end());
					}
				}
			}
		}

		Esutils::sortAndRemoveDuplicity(sHaloTarget[t]);
	}
	Esutils::mergeThreadedUniqueData(sHaloTarget);

	if (!Communication::exchangeUnknownSize(sHaloTarget[0], haloTarget, mesh._neighbours)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange elements new MPI processes.";
	}
	Esutils::mergeThreadedUniqueData(haloTarget);

	// thread x neighbor x data(oldID, new MPI process)
	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > newIDTargetSBuffet(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(targets.size()));
	std::vector<std::vector<std::pair<esglobal, esglobal> > > newIDTarget(targets.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto enodes = mesh._elems->nodes->cbegin(t);
		esglobal target;

		for (auto e = mesh._elems->distribution[t]; e < mesh._elems->distribution[t + 1]; ++e, ++enodes) {
			if (partition[e] != environment->MPIrank) {
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					auto nelems = mesh._nodes->elems->cbegin() + *n;
					for (auto ne = nelems->begin(); ne != nelems->end(); ++ne) {
						if (oldIDBoundaries[environment->MPIrank] <= *ne && *ne < oldIDBoundaries[environment->MPIrank + 1]) {
							target = partition[*ne - oldIDBoundaries[environment->MPIrank]];
						} else {
							std::pair<esglobal, esglobal> tmp(*ne, 0);
							target = std::lower_bound(haloTarget[0].begin(), haloTarget[0].end(), tmp)->second;
						}
						if (target != partition[e]) {
							newIDTargetSBuffet[t][t2i(partition[e])].push_back(std::make_pair(*ne, target));
						}
					}
				}
			}
		}
		Esutils::sortAndRemoveDuplicity(newIDTargetSBuffet[t]);
	}

	Esutils::mergeThreadedUniqueData(newIDTargetSBuffet);

	if (!Communication::sendVariousTargets(newIDTargetSBuffet[0], newIDTarget, targets)) {
		ESINFO(ERROR) << "ESPRESO internal error: exchange elements new IDs.";
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
		neighbors.push_back(begin->second);
		auto end = std::lower_bound(begin, newIDTarget[0].end(), std::make_pair(0, neighbors.back() + 1), [] (std::pair<esglobal, esglobal> &p1, const std::pair<esglobal, esglobal> &p2) { return p1.second < p2.second; });
		newID.push_back(std::move(std::vector<std::pair<esglobal, esglobal> >(begin, end)));
		begin = end;
	}

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
			newIDRequest[n][i].second = newIDBoundaries[environment->MPIrank] + it - permutation.begin();
		}
	}

//	std::vector<int> neighbors;
//	{
//		std::vector<std::vector<int> > rNeighbors(rIDsources.size());
//		#pragma omp parallel for
//		for (size_t n = 0; n < rIDsources.size(); n++) {
//			for (size_t i = 0; i < rIDsources[n].size(); i++) {
//				auto it = std::lower_bound(rNeighbors[n].begin(), rNeighbors[n].end(), rIDsources[n][i].second);
//				if (it == rNeighbors[n].end() || *it != rIDsources[n][i].second) {
//					rNeighbors[n].insert(it, rIDsources[n][i].second);
//				}
//			}
//		}
//		for (size_t n = 0; n < rIDsources.size(); n++) {
//			neighbors.insert(neighbors.end(), rNeighbors[n].begin(), rNeighbors[n].end());
//		}
//	}
//
//	std::vector<eslocal> permutation(elements->size);
//	std::iota(permutation.begin(), permutation.end(), 0);
//	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return elements->IDs->datatarray().data()[i] < elements->IDs->datatarray().data()[j]; });
//
//	#pragma omp parallel for
//	for (size_t n = 0; n < rIDsources.size(); ++n) {
//		for (size_t i = 0; i < rIDsources[n].size(); ++i) {
//			auto it = std::lower_bound(permutation.begin(), permutation.end(), rIDsources[n][i].first, [&] (esglobal i, esglobal val) {
//				return elements->IDs->datatarray().data()[i] < val;
//			});
//			rIDsources[n][i].second = newIDBoundaries[environment->MPIrank] + it - permutation.begin();
//		}
//	}
//
//	std::vector<int> nnn;
//	if (!Communication::exchangeKnownSize(rIDsources, sIDsources[0], nnn)) {
//		ESINFO(ERROR) << "ESPRESO internal error: get new elements IDs.";
//	}

	// printVector(partition);

	// printElement(elements);
	// printNodes(nodes);

	for (int rank = 0; rank < environment->MPIsize; rank++) {
		if (rank == environment->MPIrank) {
			std::cout << "RANK: " << rank << "\n";
			for (size_t t = 0; t < elements->distribution.size() - 1; t++) {
				auto ID = elements->IDs->datatarray().data();
				for (size_t e = elements->distribution[t]; e < elements->distribution[t + 1]; ++e) {
					std::cout << ID[e] << " -> " << newIDBoundaries[environment->MPIrank] + e << "\n";
				}
			}
			for (size_t n = 0; n < newIDTarget.size(); n++) {
				for (size_t i = 0; i < newIDTarget[n].size(); i++) {
					std::cout << newIDTarget[n][i].first << " -> ??\n";// << sNewID[0][n][i].second << "\n";
				}
			}
		}
		MPI_Barrier(environment->MPICommunicator);
	}
	MPI_Barrier(environment->MPICommunicator);

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "Transformation::exchanging elements finished.";
}

