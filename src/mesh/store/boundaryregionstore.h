
#ifndef SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_
#define SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_

#include <vector>
#include <string>

#include "../intervals/regioninterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct BoundaryRegionStore {

	std::string name;

	serializededata<eslocal, eslocal>* faces;
	serializededata<eslocal, eslocal>* edges;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, Element*>* facepointers;
	serializededata<eslocal, Element*>* edgepointers;

	std::vector<RegionInterval> facesIntervals;
	std::vector<RegionInterval> edgesIntervals;
	std::vector<RegionInterval> nodesIntervals;

	BoundaryRegionStore(const std::string &name);
	~BoundaryRegionStore();
};

}


#endif /* SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_ */
