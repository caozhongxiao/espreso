
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
			esint indexBase,
			size_t domains,
			const std::vector<esint> &eType,
			std::vector<std::vector<esint> > &eNodes,
			std::vector<std::vector<esint> > &eDOFs,
			std::vector<std::vector<double> > &eMatrices,
			esint dirichletSize,
			esint *dirichletIndices,
			double *dirichletValues,
			std::vector<int> &neighbours,
			size_t size, const esint *l2g) {}

protected:
	API(Mesh &mesh, esint offset, size_t domains): _mesh(mesh), _offset(offset), _domains(domains) {}

	void points(const std::vector<std::vector<esint> > &eNodes, size_t DOFsSize) {}
	void elements(const std::vector<esint> &eType, std::vector<std::vector<esint> > &eNodes, std::vector<std::vector<esint> > &eDOFs, std::vector<std::vector<double> > &eMatrices) {}
	void dirichlet(size_t dirichletSize, esint *dirichletIndices, double *dirichletValues) {}
	void clusterBoundaries(std::vector<int> &neighbours, size_t size, const esint *l2g) {}

	Mesh &_mesh;
	esint _offset;
	size_t _domains;
};

}
}




#endif /* INPUT_API_API_H_ */
