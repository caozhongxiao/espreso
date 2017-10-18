
#ifndef SRC_NEWMESH_NEWMESH_H_
#define SRC_NEWMESH_NEWMESH_H_

#include <vector>

namespace espreso {

class Mesh;
struct ElementStore;
class NewElement;
template <typename TEBoundaries, typename TEData> class serializededata;
template <typename TEBoundaries, typename TEData> struct serializededatainterval;

class MeshDomain {

public:
	int cluster;
	std::vector<int> neighbours;

	serializededatainterval<eslocal, eslocal> *elements;
};

class NewMesh {

	friend class Transformation;
public:
	NewMesh(Mesh &mesh);

	serializededata<eslocal, MeshDomain*> *_domains;

// protected:
	ElementStore *_nodes, *_edges, *_faces, *_elems, *_halo;

	// TODO: ElementStore -> RegionStore
	ElementStore *_processesCommonBoundary;

	std::vector<int> _neighbours;

	std::vector<NewElement*> _eclasses;
};

}



#endif /* SRC_NEWMESH_NEWMESH_H_ */
