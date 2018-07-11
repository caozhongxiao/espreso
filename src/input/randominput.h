
#ifndef SRC_INPUT_RANDOMINPUT_H_
#define SRC_INPUT_RANDOMINPUT_H_

#include "input.h"
#include "sfc/hilbertcurve.h"

// if SFCDEPTH > 10 => set buckets to size_t
#define SFCDEPTH 10

namespace espreso {

class RandomInput: public Input {

public:
	static void buildMesh(const ECFRoot &configuration, PlainMeshData &meshData, Mesh &mesh);

protected:
	RandomInput(const ECFRoot &configuration, PlainMeshData &dMesh, Mesh &mesh);

	void assignNBuckets();
	void assignEBuckets();

	void clusterize();
	void linkup();

	void polish();

	HilbertCurve _sfc;

	std::vector<eslocal> _nIDs;
	std::vector<uint> _nBuckets, _eBuckets;

	std::vector<eslocal> _eregions, _nregions;
	size_t _eregsize, _nregsize;

	// distribution across processes
	std::vector<uint> _bucketsBorders;
};

}


#endif /* SRC_INPUT_RANDOMINPUT_H_ */
