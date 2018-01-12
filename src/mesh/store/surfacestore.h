
#ifndef SRC_MESH_STORE_SURFACESTORE_H_
#define SRC_MESH_STORE_SURFACESTORE_H_

#include <cstddef>
#include <vector>

#include "../../basis/containers/point.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct SurfaceStore {

	serializededata<eslocal, eslocal>* triangles;
	serializededata<eslocal, eslocal>* elements;

	serializededata<eslocal, Point>* coordinates;

	std::vector<eslocal> tdistribution, edistribution, cdistribution;

	serializededata<eslocal, Element*>* epointers;
	std::vector<eslocal> ecounters;

	SurfaceStore();
	~SurfaceStore();
};

}


#endif /* SRC_MESH_STORE_SURFACESTORE_H_ */
