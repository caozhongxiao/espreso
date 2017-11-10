
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include <cstring>
#include <vector>

namespace espreso {

struct ECFConfiguration;
struct MaterialConfiguration;
class OldMesh;
class Coordinates;
class OldElement;
class Evaluator;
class Region;

namespace input {

class OldLoader {

public:
	static void load(const ECFConfiguration &configuration, OldMesh &mesh, size_t index, size_t size);

	void fill();

	virtual void points(Coordinates &coordinates) = 0;
	virtual void elements(std::vector<size_t> &bodies, std::vector<OldElement*> &elements, std::vector<OldElement*> &faces, std::vector<OldElement*> &edges) = 0;
	virtual void materials(std::vector<MaterialConfiguration*> &materials) {};
	virtual void neighbours(std::vector<OldElement*> &nodes, std::vector<int> &neighbours, const std::vector<OldElement*> &faces, const std::vector<OldElement*> &edges) = 0;
	virtual void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region*> &regions,
			std::vector<OldElement*> &elements,
			std::vector<OldElement*> &faces,
			std::vector<OldElement*> &edges,
			std::vector<OldElement*> &nodes) = 0;

	virtual void open() {};
	virtual void close() {};

	virtual bool partitiate(const std::vector<OldElement*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<OldElement*> > &fixPoints, std::vector<OldElement*> &corners) = 0;

protected:
	OldLoader(OldMesh &mesh): mesh(mesh) {}
	virtual ~OldLoader() {};

	OldMesh &mesh;
};

}
}



#endif /* INPUT_LOADER_H_ */
