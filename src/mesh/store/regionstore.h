
#ifndef SRC_MESH_STORE_REGIONSTORE_H_
#define SRC_MESH_STORE_REGIONSTORE_H_

#include <cstddef>
#include <vector>
#include <string>

#include "einterval.h"
#include "../transformation/tflags.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct RegionStore {

	std::string name;
	TFlags::ELEVEL etype;

	serializededata<eslocal, eslocal>* elems;
	serializededata<eslocal, eslocal>* faces;
	serializededata<eslocal, eslocal>* edges;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, Element*>* facepointers;
	serializededata<eslocal, Element*>* edgepointers;

	std::vector<EInterval> elemsIntervals;
	std::vector<EInterval> facesIntervals;
	std::vector<EInterval> edgesIntervals;
	std::vector<EInterval> nodesIntervals;

	std::vector<std::vector<EInterval> > domainElemsIntervals;
	std::vector<std::vector<EInterval> > domainFacesIntervals;
	std::vector<std::vector<EInterval> > domainEdgesIntervals;
	std::vector<std::vector<EInterval> > domainNodesIntervals;

	RegionStore(const std::string &name, TFlags::ELEVEL etype);
	~RegionStore();
};

}



#endif /* SRC_MESH_STORE_REGIONSTORE_H_ */
