
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include "mpi.h"

#include <string>
#include <vector>

namespace espreso {

struct ECFRoot;
struct OutputConfiguration;
struct MaterialConfiguration;

struct Statistics;
struct ElementStore;
struct ElementData;
struct NodeStore;
struct NodeData;

struct ElementsRegionStore;
struct BoundaryRegionStore;
struct ElementsRegionsIntersectionStore;
struct BoundaryRegionsIntersectionStore;
struct FETIDataStore;
struct SurfaceStore;
struct ContactStore;

class MeshPreprocessing;
class Element;

class ResultStore;

class Mesh {

	friend class MeshPreprocessing;
public:
	static Element* edata;
	static void init();
	static void destroy();

	Mesh();
	~Mesh();
	void update();
	void printMeshStatistics();
	void printDecompositionStatistics();

	ElementsRegionStore* allElements()
	{
		return eregion("ALL_ELEMENTS");
	}

	BoundaryRegionStore* allNodes()
	{
		return bregion("ALL_NODES");
	}

	ElementsRegionStore* eregion(const std::string &name);
	ElementsRegionsIntersectionStore* ieregion(const std::string &name);
	BoundaryRegionStore* bregion(const std::string &name);
	BoundaryRegionsIntersectionStore* ibregion(const std::string &name);

	bool onAllElements(const std::string &eregion) const;

	bool hasPhaseChange() const;

	void storeMesh();
	void storeSolution();

	size_t dimension;
	size_t preferedDomains;
	size_t uniformDecomposition;

	ElementStore* elements;
	NodeStore* nodes;

	std::vector<ElementsRegionStore*> elementsRegions;
	std::vector<BoundaryRegionStore*> boundaryRegions;

	std::vector<ElementsRegionsIntersectionStore*> elementsRegionsIntersections;
	std::vector<BoundaryRegionsIntersectionStore*> boundaryRegionsIntersections;

	FETIDataStore *FETIData;

	ElementStore *halo;

	SurfaceStore *surface;
	SurfaceStore *domainsSurface;

	ContactStore *contacts;

	MeshPreprocessing *preprocessing;

	std::vector<int> neighbours;
	std::vector<int> neighboursWithMe;

	std::vector<const MaterialConfiguration*> materials;

	ResultStore *store;

	bool _withGUI;
};

}



#endif /* SRC_MESH_MESH_H_ */
