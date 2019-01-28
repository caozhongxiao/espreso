
#ifndef SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_
#define SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_

#include <cstddef>
#include <vector>
#include <string>

#include "mesh/intervals/elementsinterval.h"
#include "mesh/intervals/processinterval.h"

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct Element;

struct BoundaryRegionStore {

	std::string name;

	std::vector<size_t> distribution;

	int dimension;
	double area;

	esint uniqueOffset;
	esint uniqueSize;
	esint uniqueTotalSize;

	serializededata<esint, esint>* procNodes;
	serializededata<esint, esint>* triangles;
	serializededata<esint, esint>* nodes;

	serializededata<esint, Element*>* epointers;

	std::vector<ElementsInterval> eintervals;
	std::vector<esint> eintervalsDistribution;
	std::vector<ProcessInterval> nintervals;

	std::vector<esint> ecounters;

	size_t packedSize() const;
	void pack(char* &p) const;
	void unpack(const char* &p);

	void permute(const std::vector<esint> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution);

	BoundaryRegionStore(const std::string &name, std::vector<Element*> &eclasses);
	~BoundaryRegionStore();

private:
	std::vector<Element*> &_eclasses;
};

struct BoundaryRegionsIntersectionStore: public BoundaryRegionStore {

	std::vector<BoundaryRegionStore*> regions;

	BoundaryRegionsIntersectionStore(const std::string &name, std::vector<Element*> &eclasses)
	: BoundaryRegionStore(name, eclasses) {}
};

}


#endif /* SRC_MESH_STORE_BOUNDARYREGIONSTORE_H_ */
