
#ifndef SRC_INPUT_RANDOMINPUT_H_
#define SRC_INPUT_RANDOMINPUT_H_

#include "input.h"
#include "sfc/hilbertcurve.h"

// if SFCDEPTH > 10 => set buckets to size_t
#define SFCDEPTH 10

namespace espreso {

class RandomInput: public Input {

public:
	static void buildMesh(PlainMeshData &meshData, Mesh &mesh);

protected:
	RandomInput(PlainMeshData &dMesh, Mesh &mesh);

	void assignNBuckets();
	void assignEBuckets();

	void clusterize();
	void linkup();
	void exchangeBoundary();

	void polish();

	HilbertCurve _sfc;

	std::vector<esint> _nIDs;
	std::vector<uint> _nBuckets, _eBuckets;

	// distribution across processes
	std::vector<uint> _bucketsBorders;
	std::vector<int> _sfcNeighbors;
};

}


#endif /* SRC_INPUT_RANDOMINPUT_H_ */
