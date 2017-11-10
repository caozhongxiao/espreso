
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include <vector>

namespace espreso {

class OldMesh;
struct DomainStore;
struct BoundaryStore;
struct ElementStore;
struct RegionStore;
class Element;

class Mesh {

	friend class Transformation;
public:
	Mesh(OldMesh &mesh);
	void load();

// protected:
	ElementStore *_nodes, *_edges, *_faces, *_elems, *_halo;

	DomainStore *_domains;
	BoundaryStore *_domainsBoundaries;
	BoundaryStore *_processBoundaries;

	std::vector<RegionStore*> _regions;

	std::vector<int> _neighbours;

	std::vector<Element*> _eclasses;

	OldMesh &mesh;
};

}



#endif /* SRC_MESH_MESH_H_ */
