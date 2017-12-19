
#ifndef SRC_MESH_STORE_SHAREDINTERFACESTORE_H_
#define SRC_MESH_STORE_SHAREDINTERFACESTORE_H_

#include <cstddef>
#include <utility>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct SharedInterfaceStore {

	serializededata<eslocal, eslocal>* nodes;

	std::vector<eslocal> nodeDistribution;
	std::vector<std::pair<eslocal, eslocal> > domains;

	SharedInterfaceStore();
	~SharedInterfaceStore();
};

}



#endif /* SRC_MESH_STORE_SHAREDINTERFACESTORE_H_ */
