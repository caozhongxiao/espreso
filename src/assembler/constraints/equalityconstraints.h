
#ifndef SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_

#include <utility>
#include <vector>
#include <map>
#include <cstddef>
#include <functional>

namespace espreso {

class Instance;
class Mesh;
struct ECFExpression;

class Region;
struct Step;
class SparseMatrix;

struct EqualityConstraints
{
	EqualityConstraints(Instance &instance, Mesh &mesh, const std::map<std::string, ECFExpression> &dirichlet, size_t DOFs, bool withRedundantMultiplier, bool withScaling);
	void update(const std::map<std::string, ECFExpression> &dirichlet, size_t DOFs, bool withRedundantMultiplier, bool withScaling);

	void B1DirichletInsert(const Step &step);
	void B1GlueElements(const Step &step);

	void B0Kernels(const std::vector<SparseMatrix> &kernels);
	void B0Corners();

	void insertMortarGluingToB1(const Step &step, const std::string &master, const std::string &slave);

protected:
	eslocal computeIntervalsOffsets(std::function<eslocal(eslocal)> getsize, std::function<void(eslocal, eslocal)> setsize);

	Instance &_instance;
	Mesh &_mesh;
	size_t _DOFs;
	esglobal _dirichletSize, _gluingSize;
	bool _withRedundantMultipliers, _withScaling;

	std::vector<eslocal> _mergedDirichletOffset;
	std::vector<std::vector<eslocal> > _mergedDirichletIndices;
	std::vector<std::vector<double> > _mergedDirichletValues;

	std::vector<eslocal> _intervalGluingOffset;
};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
