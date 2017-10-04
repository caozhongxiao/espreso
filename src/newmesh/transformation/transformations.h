
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

	static void computeElementCenters(NewMesh &mesh);

	static void reclusterize(NewMesh &mesh);


	static void partitiate(NewMesh &mesh, esglobal parts, TFlags::SEPARATE separate);
private:
	static size_t level;

	static void _distributeDualGraph(NewMesh &mesh, std::vector<esglobal> &edistribution, ElementStore *store, esglobal *partition, esglobal *permutation);
	static bool _checkContinuity(NewMesh &mesh, ElementStore *store);

	static void _tryrepartition(NewMesh &mesh, esglobal *permutation);
	static void _distributeNewMesh(NewMesh &mesh, ElementStore *store, esglobal *permutation);
	static void _permuteNewMesh(NewMesh &mesh, esglobal *permutation);

};

}



#endif /* SRC_MESH_TRANSFORMATION_TRANSFORMATIONS_H_ */
