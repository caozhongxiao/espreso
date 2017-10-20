
#include "transformations.h"

#include "../newmesh.h"
#include "../elements/elementstore.h"

#include "../../basis/point/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/logging/logging.h"

#include "../../config/ecf/environment.h"

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



