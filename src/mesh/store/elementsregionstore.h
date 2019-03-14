
#ifndef SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_
#define SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_

#include <vector>
#include <string>

#include "mesh/intervals/elementsinterval.h"
#include "mesh/intervals/processinterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct SurfaceStore;

struct ElementsRegionStore {

	std::string name;

	serializededata<esint, esint>* elements;
	serializededata<esint, esint>* uniqueElements;

	serializededata<esint, esint>* nodes;

	std::vector<ElementsInterval> eintervals;
	std::vector<ElementsInterval> ueintervals;
	std::vector<ProcessInterval> nintervals;

	esint uniqueOffset;
	esint uniqueSize;
	esint uniqueTotalSize;

	std::vector<esint> ecounters;

	SurfaceStore *surface;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	ElementsRegionStore(const std::string &name);
	~ElementsRegionStore();
};

}


#endif /* SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_ */
