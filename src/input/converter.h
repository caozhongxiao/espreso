
#ifndef SRC_INPUT_CONVERTER_H_
#define SRC_INPUT_CONVERTER_H_

#include "sfc/hilbertcurve.h"
#include "../basis/containers/point.h"

#include <cstddef>
#include <string>
#include <vector>

// if SFCDEPTH > 10 => set buckets to size_t
#define SFCDEPTH 9

namespace espreso {

class ECFRoot;
class Mesh;

struct EData {
	eslocal id;
	int etype;
	int body;
	int material;
};

struct MeshERegion {
	std::string name;
	std::vector<eslocal> elements;
	eslocal min, max;
};

struct MeshBRegion {
	std::string name;
	std::vector<eslocal> esize, enodes;
	std::vector<EData> edata;
	eslocal min, max;
};

struct MeshNRegion {
	std::string name;
	std::vector<eslocal> nodes;
};

struct DistributedMesh {
	std::vector<eslocal> nIDs;
	std::vector<Point> coordinates;

	std::vector<eslocal> esize, enodes;
	std::vector<EData> edata;

	std::vector<MeshERegion> eregions;
	std::vector<MeshBRegion> bregions;
	std::vector<MeshNRegion> nregions;
};

class Converter {

public:
	static void load(const ECFRoot &configuration, Mesh &mesh, int MPIrank, int MPIsize);
	static void loadDistributedMesh(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh);

protected:
	Converter(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh);

	void balance();
	void balanceNodes();
	void balancePermutedNodes();
	void balanceElements();
	void balancePermutedElements();

	void assignNBuckets();
	void assignEBuckets();

	void clusterize();
	void linkup();

	void fillElements();

	const ECFRoot &_configuration;
	DistributedMesh &_dMesh;
	Mesh &_mesh;

	HilbertCurve _sfc;

	std::vector<size_t> _nDistribution, _eDistribution;
	std::vector<eslocal> _nIDs;
	std::vector<uint> _nBuckets, _eBuckets;

	// distribution across processes
	std::vector<uint> _bucketsBorders;
};

}



#endif /* SRC_INPUT_CONVERTER_H_ */
