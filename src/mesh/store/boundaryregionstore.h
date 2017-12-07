
#ifndef SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_
#define SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_

#include <cstddef>
#include <vector>
#include <string>

#include "../intervals/elementsinterval.h"
#include "../intervals/processinterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct BoundaryRegionStore {

	std::string name;

	std::vector<size_t> distribution;

	eslocal dimension;

	eslocal uniqueOffset;
	eslocal uniqueSize;
	eslocal uniqueTotalSize;

	serializededata<eslocal, eslocal>* elements;
	serializededata<eslocal, eslocal>* nodes;

	serializededata<eslocal, Element*>* epointers;

	std::vector<ElementsInterval> eintervals;
	std::vector<ProcessInterval> nintervals;

	std::vector<eslocal> ecounters;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	void permute(const std::vector<eslocal> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);

	BoundaryRegionStore(const std::string &name, std::vector<Element*> &eclasses);
	~BoundaryRegionStore();

private:
	std::vector<Element*> &_eclasses;
};

}


#endif /* SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_ */
