
#ifndef SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_
#define SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_

#include <cstddef>
#include <vector>
#include <string>

#include "../intervals/processinterval.h"
#include "../intervals/domaininterval.h"

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

	std::vector<ProcessInterval> facesIntervals;
	std::vector<ProcessInterval> edgesIntervals;
	std::vector<ProcessInterval> nodesIntervals;

	std::vector<std::vector<DomainInterval> > domainFacesIntervals;
	std::vector<std::vector<DomainInterval> > domainEdgesIntervals;
	std::vector<std::vector<DomainInterval> > domainNodesIntervals;

	BoundaryRegionStore(const std::string &name);
	~BoundaryRegionStore();
};

}


#endif /* SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_ */
