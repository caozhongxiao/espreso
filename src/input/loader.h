
#ifndef SRC_INPUT_LOADER_H_
#define SRC_INPUT_LOADER_H_

#include "../basis/containers/point.h"

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

class ECFConfiguration;
class Mesh;

struct EData {
	int etype;
	int body;
	int material;
};

struct MeshERegion {
	std::string name;
	std::vector<eslocal> elements;
};

struct MeshBRegion {
	std::string name;
	std::vector<eslocal> esize, enodes, etypes;
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
	static void load(const ECFConfiguration &configuration, Mesh &mesh, int MPIrank, int MPIsize);
	static void loadDistributedMesh(DistributedMesh &dMesh, Mesh &mesh, bool shrinkIndices);

protected:
	Loader(DistributedMesh &dMesh, Mesh &mesh, bool shrinkIndices);

	void distributeMeshWithShrinking();
	void distributeMesh();

	void fillElements();
	void fillCoordinates();
	void addNodeRegions();
	void addBoundaryRegions();

	void shrinkIndices();

	DistributedMesh &_dMesh;
	Mesh &_mesh;

	std::vector<size_t> _nDistribution, _eDistribution;
	std::vector<int> _targetRanks;
	std::vector<std::vector<eslocal> > _rankNodeMap;
};

}

#endif /* SRC_INPUT_LOADER_H_ */
