
#ifndef SRC_OUTPUT_MESHINFO_H_
#define SRC_OUTPUT_MESHINFO_H_

#include <vector>
#include <map>
#include <string>

#include "regiondata.h"

namespace espreso {

class Mesh;
struct Region;
struct Solution;
struct Point;

namespace output {

class MeshInfo {

	friend class ResultStore;

public:
	MeshInfo(const Mesh *mesh): _mesh(mesh), _body(-1), _region(NULL) {};
	MeshInfo(const Mesh *mesh, size_t body): _mesh(mesh), _body(body), _region(NULL) {};
	MeshInfo(const Mesh *mesh, const Region *region): _mesh(mesh), _body(-1), _region(region) {};

	size_t regions() const { return _regions.size(); }
	const RegionData& region(size_t r) const { return _regions[r]; }

	virtual MeshInfo* deriveRegion(const Region *region) const =0;
	virtual MeshInfo* copyWithoutMesh() const =0;

	virtual ~MeshInfo() {};

	virtual void addSettings(size_t step) =0;
	virtual void addSolution(const std::vector<Solution*> &solution) =0;
	virtual void addGeneralInfo() = 0;

	void clearData();

	virtual bool isShrunk() const =0;
	virtual bool distributed() const =0;
	virtual Point shrink(const Point &p, eslocal domain) const =0;

protected:
	const Mesh *_mesh;

	std::vector<RegionData> _regions;
	size_t _body;
	const Region *_region;
};

}
}




#endif /* SRC_OUTPUT_MESHINFO_H_ */