
#ifndef INPUT_API_API_H_
#define INPUT_API_API_H_

#include "../loader.h"

namespace espreso {
namespace input {

class API: public Loader {

public:
	static void load(Mesh &mesh, const std::vector<std::vector<eslocal> > &eIndices, std::vector<eslocal> &neighbours, size_t size, const esglobal *ids)
	{
		ESINFO(OVERVIEW) << "Set mesh through API";
		API api(mesh, eIndices, neighbours, size, ids);
		api.fill();
	}

protected:
	// TODO: elements with various DOFS
	API(Mesh &mesh, const std::vector<std::vector<eslocal> > &eIndices, std::vector<eslocal> &neighbours, size_t size, const esglobal *ids)
	: Loader(mesh), _DOFs(3), _eIndices(eIndices), _neighbours(neighbours), _size(size), _ids(ids) { };

	void points(Coordinates &coordinates, size_t &DOFs);
	void elements(std::vector<Element*> &elements);
	void materials(std::vector<Material> &materials) { }; // unimportant
	void settings(std::vector<Evaluator*> &evaluators, std::vector<Element*> &elements, Coordinates &coordinates) {};
	void clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours);

	void fixPoints(std::vector<std::vector<eslocal> > &fixPoints) { }

private:
	size_t _DOFs;
	const std::vector<std::vector<eslocal> > &_eIndices;
	std::vector<eslocal> &_neighbours;
	size_t _size;
	const esglobal *_ids;
};

}
}




#endif /* INPUT_API_API_H_ */
