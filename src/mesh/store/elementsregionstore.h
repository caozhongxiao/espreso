
#ifndef SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_
#define SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_

#include <vector>
#include <string>

#include "../intervals/regioninterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct ElementsRegionStore {

	std::string name;

	serializededata<eslocal, eslocal>* elements;
	std::vector<RegionInterval> intervals;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	ElementsRegionStore(const std::string &name);
	~ElementsRegionStore();
};

}


#endif /* SRC_MESH_STORE_ELEMENTSREGIONSTORE_H_ */
