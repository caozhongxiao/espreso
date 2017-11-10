
#ifndef SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_

#include <vector>
#include <cstddef>
#include <functional>

namespace espreso {

class OldElement;
class Region;
class OldMesh;
enum class Property;
class Instance;
struct Step;
class SparseMatrix;

struct EqualityConstraints
{
	EqualityConstraints(
			Instance &instance,
			const OldMesh &mesh,
			const std::vector<OldElement*> &gluedElements,
			const std::vector<OldElement*> &gluedInterfaceElements,
			const std::vector<Property> &gluedDOFs,
			const std::vector<size_t> &gluedDOFsMeshOffsets,
			bool interfaceElementContainsGluedDOFs = false,
			bool dirichletSetByArrayEvaluator = false);

	void insertDirichletToB1(const Step &step, bool withRedundantMultiplier);
	void updateDirichletValuesInB1(const Step &step, bool withRedundantMultiplier);
	void insertElementGluingToB1(const Step &step, bool withRedundantMultiplier, bool withScaling);

	void insertMortarGluingToB1(const Step &step, const std::string &master, const std::string &slave);

	void insertCornersGluingToB0();
	void insertKernelsGluingToB0(const std::vector<SparseMatrix> &kernels);

protected:
	void goThroughDirichlet(
			size_t threads, const std::vector<size_t> &distribution,
			const Step &step, bool withRedundantMultiplier,
			std::function<void(size_t thread, eslocal domain, eslocal DOF, double value)> fnc);
	std::vector<esglobal> computeLambdasID(const Step &step, bool withRedundantMultiplier);

	Instance &_instance;
	const OldMesh &_mesh;
	const std::vector<OldElement*> &_gluedElements;
	const std::vector<OldElement*> &_gluedInterfaceElements;
	const std::vector<Property> &_gluedDOFs;
	const std::vector<size_t> &_gluedDOFsMeshOffsets;
	bool _interfaceElementContainsGluedDOFs;
	bool _dirichletSetByArrayEvaluator;
};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
