
#ifndef SRC_NEWMESH_NEWMESH_H_
#define SRC_NEWMESH_NEWMESH_H_

namespace espreso {

class Mesh;

struct ElementStore;

class NewMesh {

public:
	NewMesh(Mesh &mesh);

protected:
	ElementStore* _elements;
};

}



#endif /* SRC_NEWMESH_NEWMESH_H_ */
