
#ifndef SRC_MESH_STORE_SURFACESTORE_H_
#define SRC_MESH_STORE_SURFACESTORE_H_

#include <cstddef>
#include <vector>

#include "../../basis/containers/point.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct SurfaceStore {

	serializededata<eslocal, eslocal>* triangles;
	serializededata<eslocal, Point>* coordinates;

	std::vector<eslocal> edistribution, cdistribution;

	SurfaceStore();
	~SurfaceStore();
};

}


#endif /* SRC_MESH_STORE_SURFACESTORE_H_ */
