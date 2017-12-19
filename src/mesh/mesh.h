
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include "mpi.h"

#include <string>
#include <vector>

namespace espreso {

struct ECFConfiguration;
struct OutputConfiguration;
struct MaterialConfiguration;
struct Step;

struct Statistics;
struct ElementStore;
struct ElementData;
struct NodeStore;
struct NodeData;

struct ElementsRegionStore;
struct BoundaryRegionStore;
struct SharedInterfaceStore;

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

	void computeNodeStatistic(const NodeData *data, const ElementsRegionStore* region, Statistics *statistics, MPI_Comm communicator) const;
	void computeNodeStatistic(const NodeData *data, const BoundaryRegionStore* region, Statistics *statistics, MPI_Comm communicator) const;
	void computeGatheredNodeStatistic(const NodeData *data, const ElementsRegionStore* region, Statistics *statistics, MPI_Comm communicator) const;
	void computeGatheredNodeStatistic(const NodeData *data, const BoundaryRegionStore* region, Statistics *statistics, MPI_Comm communicator) const;

	void computeElementStatistic(const ElementData *data, const ElementsRegionStore* region, Statistics *statistics, MPI_Comm communicator) const;

	ElementsRegionStore* eregion(const std::string &name);
	BoundaryRegionStore* bregion(const std::string &name);

	ElementStore* elements;
	NodeStore* nodes;

	std::vector<ElementsRegionStore*> elementsRegions;
	std::vector<BoundaryRegionStore*> boundaryRegions;

	SharedInterfaceStore *sharedInterface;

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
