
#ifndef SRC_NEWMESH_STORE_BOUNDARYSTORE_H_
#define SRC_NEWMESH_STORE_BOUNDARYSTORE_H_

#include <cstddef>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct NewElement;

struct BoundaryStore {

	serializededata<eslocal, eslocal>* clusterfaces;
	serializededata<eslocal, eslocal>* localfaces;
	serializededata<eslocal, eslocal>* edges;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, NewElement*>* facepointers;
	serializededata<eslocal, NewElement*>* edgepointers;

	BoundaryStore();
	~BoundaryStore();
};

}


#endif /* SRC_NEWMESH_STORE_BOUNDARYSTORE_H_ */
