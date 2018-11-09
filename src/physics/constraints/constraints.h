
#ifndef SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_

#include <utility>
#include <vector>
#include <map>
#include <cstddef>
#include <functional>

namespace espreso {

class Instance;
class Mesh;
template <typename TValue> struct RegionMap;
struct ECFExpression;
struct ECFExpressionVector;
struct ECFExpressionOptionalVector;

class Region;
struct Step;
class SparseMatrix;
struct BoundaryRegionStore;
class Evaluator;

struct Constraints
{
	Constraints(Instance &instance, Mesh &mesh, const RegionMap<ECFExpression> &dirichlet, bool withRedundantMultiplier, bool withScaling);
	Constraints(Instance &instance, Mesh &mesh, const RegionMap<ECFExpressionOptionalVector> &dirichlet, int DOFs, bool withRedundantMultiplier, bool withScaling);
	void update(const RegionMap<ECFExpression> &dirichlet, bool withRedundantMultiplier, bool withScaling);
	void update(const RegionMap<ECFExpressionOptionalVector> &dirichlet, bool withRedundantMultiplier, bool withScaling);

	void B1DirichletInsert(const Step &step);
	void B1DirichletUpdate(const Step &step);
	void B1GlueElements();
	void B1DuplicityUpdate();

	void B0Kernels(const std::vector<SparseMatrix> &kernels);
	void B0Corners();

	void B1ContactInsert(const Step &step, BoundaryRegionStore *region, const ECFExpressionVector &normals, const ECFExpressionVector &gap);

	void insertMortarGluingToB1(const Step &step, const std::string &master, const std::string &slave);

protected:
	eslocal computeIntervalsOffsets(std::function<eslocal(eslocal)> getsize, std::function<void(eslocal, eslocal)> setsize);
	void update(bool withRedundantMultiplier, bool withScaling);

	Instance &_instance;
	Mesh &_mesh;
	esglobal _dirichletSize, _gluingSize, _contactSize;
	bool _withRedundantMultipliers, _withScaling;

	int _DOFs;

	// DOF x CONDITIONS
	std::vector<std::vector<std::pair<BoundaryRegionStore*, Evaluator*> > > _dirichlet;

	// DOF x INTERVAL x NODES
	std::vector<std::vector<std::vector<eslocal> > > _intervalDirichletNodes;

	std::vector<eslocal> _domainDirichletSize;

	std::vector<eslocal> _intervalDirichletOffset;
	std::vector<eslocal> _intervalGluingOffset;
};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_ */
