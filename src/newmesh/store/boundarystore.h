
#ifndef SRC_NEWMESH_STORE_BOUNDARYSTORE_H_
#define SRC_NEWMESH_STORE_BOUNDARYSTORE_H_

#include <cstddef>
#include <vector>

#include "einterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct NewElement;

struct BoundaryStore {

	serializededata<eslocal, eslocal>* elems;
	serializededata<eslocal, eslocal>* faces;
	serializededata<eslocal, eslocal>* edges;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, NewElement*>* facepointers;
	serializededata<eslocal, NewElement*>* edgepointers;

	std::vector<EInterval> elemsIntervals;
	std::vector<EInterval> facesIntervals;
	std::vector<EInterval> edgesIntervals;
	std::vector<EInterval> nodesIntervals;

	BoundaryStore();
	~BoundaryStore();
};

}


#endif /* SRC_NEWMESH_STORE_BOUNDARYSTORE_H_ */
