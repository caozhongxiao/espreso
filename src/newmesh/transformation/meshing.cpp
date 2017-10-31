
#include "transformations.h"

#include "../newmesh.h"
#include "../elements/elementstore.h"

#include "../../basis/point/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../../config/ecf/environment.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

void Transformation::computeElementCenters(NewMesh &mesh)
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

void Transformation::computeDomainsCenters(NewMesh &mesh)
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

void Transformation::reindexNodes(NewMesh &mesh)
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
		if (n < environment->MPIrank) {
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

void Transformation::arrangeNodes(NewMesh &mesh)
{
	ESINFO(TVERBOSITY) << std::string(2 * level++, ' ') << "MESH::arrange nodes started.";

	if (mesh._nodes->domains == NULL) {
		Transformation::assignDomainsToNodes(mesh);
	}

	size_t threads = environment->OMP_NUM_THREADS;

	std::vector<eslocal> permutation(mesh._nodes->size);
	std::iota(permutation.begin(), permutation.end(), 0);

	// sorted according the following criteria:
	// 1. on external boundary
	// 2. on internal boundary
	// 3. number of boundaries domains
	// 4. boundaries IDs
	// 5. ID
//	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
//		auto di = mesh._nodes->domains->cbegin() + i;
//		auto dj = mesh._nodes->domains->cbegin() + j;
//
//		if ((di->front() < 0 || di->back() < 0) && (dj->front() >= 0 || dj->back() >= 0)) {
//			// only first node is on boundary
//			return true;
//		}
//		if ((dj->front() < 0 || dj->back() < 0) && (di->front() >= 0 || di->back() >= 0)) {
//			// only second node is on boundary
//			return false;
//		}
//		if (di->size() == dj->size()) {
//			// the same number of boundaries
//			for (size_t d = 0; d < di->size(); d++) {
//				if ((*di)[d] != (*dj)[d]) {
//					return (*di)[d] < (*dj)[d];
//				}
//			}
//			return mesh._nodes->IDs->datatarray()[i] < mesh._nodes->IDs->datatarray()[j];
//		}
//		return di->size() < dj->size();
//	});

	// Communication::serialize([&] () { std::cout << permutation; });

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

	ESINFO(TVERBOSITY) << std::string(--level * 2, ' ') << "MESH::arrange nodes finished.";
}
