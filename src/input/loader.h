
#ifndef SRC_INPUT_LOADER_H_
#define SRC_INPUT_LOADER_H_

#include "../basis/containers/point.h"

#include <cstddef>
#include <string>
#include <vector>

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

class Loader {

public:
	static void load(const ECFRoot &configuration, Mesh &mesh, int MPIrank, int MPIsize);
	static void loadDistributedMesh(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh);

protected:
	Loader(const ECFRoot &configuration, DistributedMesh &dMesh, Mesh &mesh);

	void distributeMesh();

	void fillElements();
	void fillCoordinates();
	void addNodeRegions();
	void addBoundaryRegions();
	void addElementRegions();

	const ECFRoot &_configuration;
	DistributedMesh &_dMesh;
	Mesh &_mesh;

	std::vector<size_t> _nDistribution, _eDistribution;
	std::vector<int> _targetRanks;
	std::vector<std::vector<eslocal> > _rankNodeMap;
};

}

#endif /* SRC_INPUT_LOADER_H_ */
