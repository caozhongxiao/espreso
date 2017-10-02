
#ifndef SRC_NEWMESH_NEWMESH_H_
#define SRC_NEWMESH_NEWMESH_H_

namespace espreso {

class Mesh;
struct ElementStore;

class NewMesh {

public:
	NewMesh(Mesh &mesh);

protected:
	ElementStore *_nodes, *_edges, *_faces, *_elems;

};

}



#endif /* SRC_NEWMESH_NEWMESH_H_ */
