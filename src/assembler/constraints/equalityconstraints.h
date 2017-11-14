
#ifndef SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_

#include <vector>
#include <cstddef>
#include <functional>

namespace espreso {

class Instance;
class Mesh;
class RegionStore;

class OldElement;
class Region;
enum class Property;
struct Step;
class SparseMatrix;

struct EqualityConstraints
{
	EqualityConstraints(Instance &instance, Mesh &mesh, const std::vector<RegionStore*> &dirichlet, size_t DOFs, bool withRedundantMultiplier, bool withScaling);
	void update(const std::vector<RegionStore*> &dirichlet, size_t DOFs, bool withRedundantMultiplier, bool withScaling);

	void B1DirichletInsert(const Step &step);
	void B1GlueElements(const Step &step);


	void insertMortarGluingToB1(const Step &step, const std::string &master, const std::string &slave);
	void insertCornersGluingToB0();
	void insertKernelsGluingToB0(const std::vector<SparseMatrix> &kernels);

protected:
	Instance &_instance;
	Mesh &_mesh;
	std::vector<RegionStore*> _dirichlet;
	size_t _DOFs;
	esglobal _dirichletSize, _gluingSize;
	bool _withRedundantMultipliers, _withScaling;

	std::vector<std::vector<eslocal> > _mergedDirichletIndices;
	std::vector<std::vector<double> > _mergedDirichletValues;
};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
