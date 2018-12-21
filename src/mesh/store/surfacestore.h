
#ifndef SRC_MESH_STORE_SURFACESTORE_H_
#define SRC_MESH_STORE_SURFACESTORE_H_

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Point;
struct Element;

struct SurfaceStore {

	serializededata<esint, esint>* triangles;
	serializededata<esint, esint>* elements;

	serializededata<esint, Point>* coordinates;

	std::vector<size_t> tdistribution, edistribution, cdistribution;

	serializededata<esint, Element*>* epointers;
	std::vector<esint> ecounters;

	SurfaceStore();
	~SurfaceStore();
};

}


#endif /* SRC_MESH_STORE_SURFACESTORE_H_ */
