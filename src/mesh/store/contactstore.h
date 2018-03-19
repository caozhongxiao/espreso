
#ifndef SRC_MESH_STORE_CONTACTSTORE_H_
#define SRC_MESH_STORE_CONTACTSTORE_H_

#include <cstddef>
#include <vector>

#include "../../basis/containers/point.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct SurfaceStore;

struct ContactStore {

	double eps;
	size_t groupsize;

	SurfaceStore *surface;
	serializededata<eslocal, Point>* elements;

	serializededata<eslocal, eslocal>* closeElements;

	Point boundingBox[2], globalBox[2];

	std::vector<eslocal> filledCells;
	serializededata<eslocal, eslocal>* grid;

	std::vector<int> neighbors;
	std::vector<std::vector<eslocal> > nsurface;
	std::vector<serializededata<eslocal, Point>*> nelements;
	std::vector<std::vector<eslocal> > nfilled;
	std::vector<serializededata<eslocal, eslocal>*> ngrid;
	std::vector<serializededata<eslocal, eslocal>*> ncloseElements;

	ContactStore(SurfaceStore *surface);
	~ContactStore();
};

}


#endif /* SRC_MESH_STORE_CONTACTSTORE_H_ */
