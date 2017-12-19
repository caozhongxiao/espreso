
#ifndef SRC_MESH_STORE_FETIDATASTORE_H_
#define SRC_MESH_STORE_FETIDATASTORE_H_

#include <cstddef>
#include <utility>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct FETIDataStore {

	// B0 from kernels
	serializededata<eslocal, eslocal>* interfaceNodes;
	std::vector<eslocal> inodesDistribution;
	std::vector<std::pair<eslocal, eslocal> > inodesDomains;

	// B0 from corners

	std::vector<eslocal> corners;
	serializededata<eslocal, eslocal>* cornerDomains;

	FETIDataStore();
	~FETIDataStore();
};

}



#endif /* SRC_MESH_STORE_FETIDATASTORE_H_ */
