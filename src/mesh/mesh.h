
#ifndef SRC_MESH_MESH_H_
#define SRC_MESH_MESH_H_

#include <string>
#include <vector>
#include <functional>

namespace espreso {

class OldMesh;
struct DomainStore;
struct BoundaryStore;
struct ElementStore;
struct RegionStore;
class Element;
struct EInterval;

enum class ETYPE: int {
	ELEMENT,
	FACE,
	EDGE,
	NODE
};

class Mesh {

	friend class Transformation;
public:
	Mesh();
	void load();

	esglobal computeIntervalsOffsets(std::vector<EInterval> &intervals, std::function<eslocal(eslocal)> getsize, std::function<void(eslocal, esglobal)> setsize);

	RegionStore* region(const std::string &name);

// protected:
	ElementStore *_nodes, *_edges, *_faces, *_elems, *_halo;

	DomainStore *_domains;
	BoundaryStore *_domainsBoundaries;
	BoundaryStore *_processBoundaries;

	std::vector<RegionStore*> _regions;

	std::vector<int> _neighbours;

	std::vector<Element*> _eclasses;

	OldMesh *mesh;
};

}



#endif /* SRC_MESH_MESH_H_ */
