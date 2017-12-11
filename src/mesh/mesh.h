
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include <string>
#include <vector>

namespace espreso {

struct ECFConfiguration;
struct OutputConfiguration;
struct MaterialConfiguration;
struct Step;

struct Statistics;
struct ElementStore;
struct NodeStore;

struct ElementsRegionStore;
struct BoundaryRegionStore;

class MeshPreprocessing;
class Element;

class OldMesh;

class Mesh {

	friend class MeshPreprocessing;
public:
	Mesh(const ECFConfiguration &configuration);
	void load();
	void update();

	void initNodeData();
	void gatherNodeData();
	void computeNodeStatistic(const ElementsRegionStore* region, Statistics &statistics) const;
	void computeNodeStatistic(const BoundaryRegionStore* region, Statistics &statistics) const;
	void computeGatheredNodeStatistic(const ElementsRegionStore* region, Statistics &statistics) const;
	void computeGatheredNodeStatistic(const BoundaryRegionStore* region, Statistics &statistics) const;

	ElementsRegionStore* eregion(const std::string &name);
	BoundaryRegionStore* bregion(const std::string &name);

	ElementStore* elements;
	NodeStore* nodes;

	std::vector<ElementsRegionStore*> elementsRegions;
	std::vector<BoundaryRegionStore*> boundaryRegions;

	ElementStore *halo;

	MeshPreprocessing *preprocessing;

	std::vector<int> neighbours;
	std::vector<int> neighboursWithMe;

	std::vector<const MaterialConfiguration*> materials;

//protected:
	const ECFConfiguration &configuration;
	std::vector<Element*> _eclasses;
	OldMesh *mesh;

};

}



#endif /* SRC_MESH_MESH_H_ */
