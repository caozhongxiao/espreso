
#ifndef SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_
#define SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_

#include <cstddef>
#include <vector>

namespace espreso {

#define VERBOSITY(level) level == 0 ? PROGRESS1 : level == 1 ? PROGRESS2 : PROGRESS3

class Mesh;

class MeshPreprocessing {

	friend class Mesh;
public:

	void linkNodesAndElements();
	void exchangeHalo();

	void computeDual();
	void computeDecomposedDual(bool separateMaterials, bool separateEtype);

	void reclusterize();
	void partitiate(eslocal parts, bool separateMaterials, bool separateEtype);

protected:
	static size_t level;

	MeshPreprocessing(Mesh *mesh): _mesh(mesh) {}
	Mesh *_mesh;

private:
	void exchangeElements(const std::vector<eslocal> &partition);
	void permuteElements(const std::vector<eslocal> &permutation, const std::vector<size_t> &distribution);
};

}



#endif /* SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_ */
