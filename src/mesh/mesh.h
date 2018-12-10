
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
	Mesh(const ECFRoot &configuration, ResultStore *store, bool withGUI = false);
	void update();

	void initNodeData();
	void gatherNodeData();

	double sumSquares(const std::vector<std::vector<double> > &data, const BoundaryRegionStore* region);
	void computeGatheredNodeStatistic(const NodeData *data, const ElementsRegionStore* region, Statistics *statistics, MPI_Comm communicator) const;
	void computeGatheredNodeStatistic(const NodeData *data, const BoundaryRegionStore* region, Statistics *statistics, MPI_Comm communicator) const;

	void computeElementStatistic(const ElementData *data, const ElementsRegionStore* region, Statistics *statistics, MPI_Comm communicator) const;

	ElementsRegionStore* eregion(const std::string &name);
	ElementsRegionsIntersectionStore* ieregion(const std::string &name);
	BoundaryRegionStore* bregion(const std::string &name);
	BoundaryRegionsIntersectionStore* ibregion(const std::string &name);

	bool onAllElements(const std::string &eregion) const;

	bool hasPhaseChange() const;

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

//protected:
	const ECFRoot &configuration;
	std::vector<Element*> _eclasses;

	ResultStore *store;

	bool _withGUI;

private:
	void printMeshStatistics();
	void printDecompositionStatistics();
};

}



#endif /* SRC_MESH_MESH_H_ */
