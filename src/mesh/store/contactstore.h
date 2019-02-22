
#ifndef SRC_MESH_STORE_CONTACTSTORE_H_
#define SRC_MESH_STORE_CONTACTSTORE_H_

#include <cstddef>
#include <vector>

#include "basis/containers/point.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct SurfaceStore;

struct ContactStore {

	double eps;
	size_t groupsize;

	SurfaceStore *surface;
	serializededata<esint, Point>* elements;
	serializededata<esint, Point>* enormals;

	serializededata<esint, esint>* closeElements;

	Point boundingBox[2], globalBox[2];

	std::vector<esint> filledCells;
	serializededata<esint, esint>* grid;

	std::vector<int> neighbors;
	std::vector<std::vector<esint> > nsurface;
	std::vector<serializededata<esint, Point>*> nelements;
	std::vector<std::vector<esint> > nfilled;
	std::vector<serializededata<esint, esint>*> ngrid;
	std::vector<serializededata<esint, esint>*> ncloseElements;

	ContactStore(SurfaceStore *surface);
	~ContactStore();
};

}


#endif /* SRC_MESH_STORE_CONTACTSTORE_H_ */
