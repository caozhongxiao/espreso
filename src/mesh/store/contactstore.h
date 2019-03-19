
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

	// geometric
	std::vector<int> gneighbors;
	// neighbors
	std::vector<std::vector<esint> > gnsurface;
	std::vector<serializededata<esint, Point>*> gnelements;
	std::vector<std::vector<esint> > gnfilled;
	std::vector<serializededata<esint, esint>*> gngrid;
	std::vector<serializededata<esint, esint>*> gncloseElements;

	ContactStore(SurfaceStore *surface);
	~ContactStore();

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);
};

}


#endif /* SRC_MESH_STORE_CONTACTSTORE_H_ */
