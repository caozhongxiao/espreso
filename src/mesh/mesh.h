
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include <vector>
#include <functional>

namespace espreso {

struct ElementStore;
struct NodeStore;

struct ElementsRegionStore;
struct BoundaryRegionStore;

class MeshPreprocessing;

// OLD

class OldMesh;
struct DomainStore;
struct BoundaryStore;
struct RegionStore;
class Element;
struct EInterval;
class MaterialConfiguration;
struct ECFConfiguration;

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

	esglobal computeIntervalsOffsets(std::vector<EInterval> &intervals, std::function<eslocal(eslocal)> getsize, std::function<void(eslocal, esglobal)> setsize);

	// RegionStore* region(const std::string &name) {};

	ElementStore* elements;
	NodeStore* nodes;

	std::vector<ElementsRegionStore*> elementsRegions;
	std::vector<BoundaryRegionStore*> boundaryRegions;

	ElementStore *halo;

	MeshPreprocessing *preprocessing;

	std::vector<int> neighbours;

	///// <<<<< OLD >>>>>> //////

	DomainStore *_domains;
	BoundaryStore *_domainsBoundaries;
	BoundaryStore *_processBoundaries;

	std::vector<const MaterialConfiguration*> _materials;



//protected:
	std::vector<Element*> _eclasses;
	OldMesh *mesh;
};

}



#endif /* SRC_MESH_MESH_H_ */
