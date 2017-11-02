
#ifndef SRC_NEWMESH_STORE_DOMAINSTORE_H_
#define SRC_NEWMESH_STORE_DOMAINSTORE_H_

#include <cstddef>
#include <vector>

#include "einterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct DomainStore {

	size_t size;
	std::vector<size_t> domainDistribution;
	std::vector<size_t> domainElementBoundaries;
	std::vector<size_t> domainNodeBoundaries;

	std::vector<int> clusters;

	serializededata<eslocal, eslocal> *elems;
	serializededata<eslocal, eslocal> *nodes;

	// domain x intervals (in order: internal boundary, external boundary, internal nodes)
	std::vector<std::vector<EInterval> > nodesIntervals;

	std::vector<esglobal> gatherDomainDistribution();

	DomainStore();
	~DomainStore();
};

}

#endif /* SRC_NEWMESH_STORE_DOMAINSTORE_H_ */
