
#ifndef SRC_INPUT_LOADER_H_
#define SRC_INPUT_LOADER_H_

#include "../basis/containers/point.h"

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
	std::vector<eslocal> edist, enodes, etypes;
};

struct MeshNRegion {
	std::string name;
	std::vector<eslocal> nodes;
};

struct DistributedMesh {
	std::vector<eslocal> nIDs;
	std::vector<Point> coordinates;

	std::vector<eslocal> edist, enodes;
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

	void fillMesh();
	void addNodeRegions();
	void addBoundaryRegions();

	void shrinkIndices();

	DistributedMesh &_dMesh;
	Mesh &_mesh;
};

}

#endif /* SRC_INPUT_LOADER_H_ */
