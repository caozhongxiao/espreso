
#ifndef SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_
#define SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_

#include <vector>
#include <string>

#include "../intervals/processinterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct BoundaryRegionStore {

	std::string name;

	eslocal uniqueOffset;
	eslocal uniqueSize;
	eslocal uniqueTotalSize;

	serializededata<eslocal, eslocal>* faces;
	serializededata<eslocal, eslocal>* edges;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, Element*>* facepointers;
	serializededata<eslocal, Element*>* edgepointers;

	std::vector<ProcessInterval> facesIntervals;
	std::vector<ProcessInterval> edgesIntervals;
	std::vector<ProcessInterval> nodesIntervals;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	BoundaryRegionStore(const std::string &name, std::vector<Element*> &eclasses);
	~BoundaryRegionStore();

private:
	std::vector<Element*> &_eclasses;
};

}


#endif /* SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_ */
