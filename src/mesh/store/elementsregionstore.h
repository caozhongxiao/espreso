
#ifndef SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_
#define SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_

#include <vector>
#include <string>

#include "../intervals/elementsinterval.h"
#include "../intervals/processinterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct ElementsRegionStore {

	std::string name;

	serializededata<eslocal, eslocal>* elements;
	serializededata<eslocal, eslocal>* uniqueElements;

	serializededata<eslocal, eslocal>* nodes;

	std::vector<ElementsInterval> eintervals;
	std::vector<ElementsInterval> ueintervals;
	std::vector<ProcessInterval> nintervals;

	eslocal uniqueOffset;
	eslocal uniqueSize;
	eslocal uniqueTotalSize;

	std::vector<eslocal> ecounters;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	ElementsRegionStore(const std::string &name);
	~ElementsRegionStore();
};

struct ElementsRegionsIntersectionStore: public ElementsRegionStore {

	std::vector<ElementsRegionStore*> regions;

	ElementsRegionsIntersectionStore(const std::string &name): ElementsRegionStore(name) {}
};

}


#endif /* SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_ */
