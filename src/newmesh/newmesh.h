
#ifndef SRC_NEWMESH_NEWMESH_H_
#define SRC_NEWMESH_NEWMESH_H_

#include <vector>

namespace espreso {

class Mesh;
struct DomainStore;
struct BoundaryStore;
struct ElementStore;
class NewElement;

class NewMesh {

	friend class Transformation;
public:
	NewMesh(Mesh &mesh);

	DomainStore *_domains;
	BoundaryStore *_domainsBoundaries;
	BoundaryStore *_processBoundaries;

// protected:
	ElementStore *_nodes, *_edges, *_faces, *_elems, *_halo;

	std::vector<int> _neighbours;

	std::vector<NewElement*> _eclasses;
};

}



#endif /* SRC_NEWMESH_NEWMESH_H_ */
