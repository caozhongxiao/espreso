
#ifndef SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_
#define SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_

#include "tflags.h"

#include <vector>

namespace espreso {

class NewMesh;
class ElementStore;

#define TVERBOSITY PROGRESS1

struct Transformation {

	static void addLinkFromTo(NewMesh &mesh, TFlags::ELEVEL from, TFlags::ELEVEL to);

	static void exchangeHaloElements(NewMesh &mesh);

	static void computeDual(NewMesh &mesh);
	static void computeDecomposedDual(NewMesh &mesh, TFlags::SEPARATE separate);

	static void computeProcessesCommonBoundary(NewMesh &mesh);
	static void computeDomainsBoundaries(NewMesh &mesh);

	static void computeElementCenters(NewMesh &mesh);
	static void computeDomainsCenters(NewMesh &mesh);

	static void reclusterize(NewMesh &mesh);
	static void partitiate(NewMesh &mesh, esglobal parts, TFlags::SEPARATE separate);
	static void assignDomainsToNodes(NewMesh &mesh);
private:
	static size_t level;

	static void exchangeElements(NewMesh &mesh, const std::vector<esglobal> &partition);
	static void permuteElements(NewMesh &mesh, const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);
	static void permuteNodes(NewMesh &mesh, const std::vector<eslocal> &permutation);
};

}



#endif /* SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_ */
