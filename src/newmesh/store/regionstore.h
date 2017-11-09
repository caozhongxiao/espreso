
#ifndef SRC_NEWMESH_STORE_REGIONSTORE_H_
#define SRC_NEWMESH_STORE_REGIONSTORE_H_

#include <cstddef>
#include <vector>
#include <string>

#include "einterval.h"
#include "../transformation/tflags.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct NewElement;

struct RegionStore {

	std::string name;
	TFlags::ELEVEL etype;

	serializededata<eslocal, eslocal>* elems;
	serializededata<eslocal, eslocal>* faces;
	serializededata<eslocal, eslocal>* edges;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, NewElement*>* facepointers;
	serializededata<eslocal, NewElement*>* edgepointers;

	std::vector<std::vector<EInterval> > elemsIntervals;
	std::vector<std::vector<EInterval> > facesIntervals;
	std::vector<std::vector<EInterval> > edgesIntervals;
	std::vector<std::vector<EInterval> > nodesIntervals;

	RegionStore(const std::string &name, TFlags::ELEVEL etype);
	~RegionStore();
};

}



#endif /* SRC_NEWMESH_STORE_REGIONSTORE_H_ */
