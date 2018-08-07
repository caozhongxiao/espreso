
#ifndef INPUT_API_API_H_
#define INPUT_API_API_H_

// TODO: MESH

namespace espreso {

class Mesh;

namespace input {

class API {

public:
	static void load(
			Mesh &mesh,
			eslocal indexBase,
			size_t domains,
			const std::vector<eslocal> &eType,
			std::vector<std::vector<eslocal> > &eNodes,
			std::vector<std::vector<eslocal> > &eDOFs,
			std::vector<std::vector<double> > &eMatrices,
			eslocal dirichletSize,
			eslocal *dirichletIndices,
			double *dirichletValues,
			std::vector<int> &neighbours,
			size_t size, const eslocal *l2g) {}

protected:
	API(Mesh &mesh, eslocal offset, size_t domains): _mesh(mesh), _offset(offset), _domains(domains) {}

	void points(const std::vector<std::vector<eslocal> > &eNodes, size_t DOFsSize) {}
	void elements(const std::vector<eslocal> &eType, std::vector<std::vector<eslocal> > &eNodes, std::vector<std::vector<eslocal> > &eDOFs, std::vector<std::vector<double> > &eMatrices) {}
	void dirichlet(size_t dirichletSize, eslocal *dirichletIndices, double *dirichletValues) {}
	void clusterBoundaries(std::vector<int> &neighbours, size_t size, const eslocal *l2g) {}

	Mesh &_mesh;
	eslocal _offset;
	size_t _domains;
};

}
}




#endif /* INPUT_API_API_H_ */
