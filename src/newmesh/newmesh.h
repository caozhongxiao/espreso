
#ifndef SRC_NEWMESH_NEWMESH_H_
#define SRC_NEWMESH_NEWMESH_H_

#include "elements/newelement.h"

#include <vector>

namespace espreso {

class Mesh;
struct ElementStore;

class NewMesh {

	friend class Transformation;
public:
	NewMesh(Mesh &mesh);

// protected:
	ElementStore *_nodes, *_edges, *_faces, *_elems, *_halo;
	std::vector<int> _neighbours;

	std::vector<std::vector<NewElement> > _eclasses;
};

}



#endif /* SRC_NEWMESH_NEWMESH_H_ */
