
#ifndef SRC_NEWMESH_STORE_BOUNDARYSTORE_H_
#define SRC_NEWMESH_STORE_BOUNDARYSTORE_H_

#include <cstddef>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct NewElement;

struct BoundaryInterval {
	std::vector<int> neighbors;
	eslocal begin, end;

	BoundaryInterval(eslocal begin, eslocal end, std::vector<int> &&neighbors)
	: begin(begin), end(end), neighbors(std::move(neighbors)) {};
};

struct BoundaryStore {

	serializededata<eslocal, eslocal>* elems;
	serializededata<eslocal, eslocal>* faces;
	serializededata<eslocal, eslocal>* edges;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, NewElement*>* facepointers;
	serializededata<eslocal, NewElement*>* edgepointers;

	std::vector<BoundaryInterval> elemsIntervals;
	std::vector<BoundaryInterval> facesIntervals;
	std::vector<BoundaryInterval> edgesIntervals;
	std::vector<BoundaryInterval> nodesIntervals;

	BoundaryStore();
	~BoundaryStore();
};

}


#endif /* SRC_NEWMESH_STORE_BOUNDARYSTORE_H_ */
