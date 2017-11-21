
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include <string>
#include <vector>

namespace espreso {

struct ECFConfiguration;
struct MaterialConfiguration;

struct ElementStore;
struct NodeStore;

struct ElementsRegionStore;
struct BoundaryRegionStore;

class MeshPreprocessing;
class Element;

class OldMesh;

enum class ETYPE: int {
	ELEMENT,
	FACE,
	EDGE,
	NODE
};

class Mesh {

	friend class MeshPreprocessing;
public:
	Mesh();
	void load(const ECFConfiguration &configuration);
	void update(const ECFConfiguration &configuration);

	ElementsRegionStore* eregion(const std::string &name);
	BoundaryRegionStore* bregion(const std::string &name);

	ElementStore* elements;
	NodeStore* nodes;

	std::vector<ElementsRegionStore*> elementsRegions;
	std::vector<BoundaryRegionStore*> boundaryRegions;

	ElementStore *halo;

	MeshPreprocessing *preprocessing;

	std::vector<int> neighbours;

	std::vector<const MaterialConfiguration*> materials;

//protected:
	std::vector<Element*> _eclasses;
	OldMesh *mesh;
};

}



#endif /* SRC_MESH_MESH_H_ */
