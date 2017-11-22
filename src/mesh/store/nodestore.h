
#ifndef SRC_MESH_STORE_NODESTORE_H_
#define SRC_MESH_STORE_NODESTORE_H_

#include <cstddef>
#include <string>
#include <vector>

#include "../intervals/processinterval.h"
#include "../intervals/domaininterval.h"
#include "../intervals/gluinginterval.h"

namespace espreso {

struct Point;
template <typename TEBoundaries, typename TEData> class serializededata;

struct NodeStore {

	void store(const std::string &file);

	void permute(const std::vector<eslocal> &permutation) { permute(permutation, distribution); }
	void permute(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);

	std::vector<eslocal> gatherNodeDistribution();

	eslocal size;
	std::vector<size_t> distribution;

	serializededata<eslocal, eslocal>* IDs;
	serializededata<eslocal, eslocal>* elements;

	serializededata<eslocal, Point>* coordinates;
	serializededata<eslocal, int>* ranks;

	std::vector<eslocal> externalIntervals;
	std::vector<ProcessInterval> pintervals;
	std::vector<std::vector<DomainInterval> > dintervals;
	std::vector<std::vector<GluingInterval> > gintervals;
	serializededata<eslocal, eslocal>* idomains;
	serializededata<eslocal, int>* iprocesses;

	NodeStore();
	~NodeStore();
};

}


#endif /* SRC_MESH_STORE_NODESTORE_H_ */
