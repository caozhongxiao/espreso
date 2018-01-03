
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
struct BoundaryRegionStore;
class Evaluator;

struct EqualityConstraints
{
	EqualityConstraints(Instance &instance, Mesh &mesh, const std::map<std::string, ECFExpression> &dirichlet, bool withRedundantMultiplier, bool withScaling);
	void update(const std::map<std::string, ECFExpression> &dirichlet, bool withRedundantMultiplier, bool withScaling);

	void B1DirichletInsert(const Step &step);
	void B1DirichletUpdate(const Step &step);
	void B1GlueElements();
	void B1DuplicityUpdate();

	void B0Kernels(const std::vector<SparseMatrix> &kernels);
	void B0Corners();

	void insertMortarGluingToB1(const Step &step, const std::string &master, const std::string &slave);

protected:
	eslocal computeIntervalsOffsets(std::function<eslocal(eslocal)> getsize, std::function<void(eslocal, eslocal)> setsize);

	Instance &_instance;
	Mesh &_mesh;
	esglobal _dirichletSize, _gluingSize;
	bool _withRedundantMultipliers, _withScaling;

	// DOF x CONDITIONS
	std::vector<std::vector<std::pair<BoundaryRegionStore*, Evaluator*> > > _dirichlet;

	// DOF x INTERVAL x NODES
	std::vector<std::vector<std::vector<eslocal> > > _intervalDirichletNodes;

	std::vector<eslocal> _domainDirichletSize;

	std::vector<eslocal> _intervalDirichletOffset;
	std::vector<eslocal> _intervalGluingOffset;
};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
