
#include "transformations.h"

#include "../../config/ecf/environment.h"

#include "../elements/element.h"
#include "../elements/elementstore.h"

#include "../newmesh.h"

using namespace espreso;

void Transformation::computeElementCenters(NewMesh &mesh)
{
//	const parray<Element*> &elements = *mesh.__elements->elements;
//	size_t threads = environment->OMP_NUM_THREADS;
//
//	std::vector<size_t> distribution = mesh.__elements->distribution;
//	for (size_t i = 0; i < distribution.size(); i++) {
//		distribution[i] *= mesh.__elements->dimension;
//	}
//
//	mesh.__elements->coordinates = new parray<double>(distribution.size(), distribution.data());
//
//	parray<double> &centers = *mesh.__elements->coordinates;
//	size_t dimension = mesh.__elements->dimension;
//
//	#pragma omp parallel for
//	for (size_t t = 0; t < threads; t++) {
//		const crsiterator<eslocal> &enodes = mesh.__elements->nodesIndices->iterator(t);
//
//		while (enodes.next()) {
//			Point center;
//			for (auto n = enodes.begin(); n != enodes.end(); ++n) {
//				center += mesh.coordinates()[*n];
//			}
//			center /= enodes.size();
//
//			centers[dimension * enodes.index() + 0] = center.x;
//			centers[dimension * enodes.index() + 1] = center.y;
//			if (dimension == 3) {
//				centers[dimension * enodes.index() + 2] = center.z;
//			}
//		}
//	}
}



