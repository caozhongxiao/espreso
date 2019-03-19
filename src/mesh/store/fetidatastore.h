
#ifndef SRC_MESH_STORE_FETIDATASTORE_H_
#define SRC_MESH_STORE_FETIDATASTORE_H_

#include <cstddef>
#include <utility>
#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct FETIDataStore {

	// B0 from kernels
	serializededata<esint, esint>* interfaceNodes;
	std::vector<esint> inodesDistribution;
	std::vector<std::pair<esint, esint> > inodesDomains;

	// B0 from corners

	std::vector<esint> corners;
	serializededata<esint, esint>* cornerDomains;

	// Regularization from fix points
	std::vector<esint> surfaceFixPoints, sFixPointsDistribution;
	std::vector<esint> innerFixPoints, iFixPointsDistribution;

	FETIDataStore();
	~FETIDataStore();

	size_t packedFullSize() const;
	void packFull(char* &p) const;
	void unpackFull(const char* &p);
};

}



#endif /* SRC_MESH_STORE_FETIDATASTORE_H_ */
