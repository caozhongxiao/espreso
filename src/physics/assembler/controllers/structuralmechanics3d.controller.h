
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS3D_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS3D_CONTROLLER_H_

#include "structuralmechanics.controller.h"

namespace espreso {

struct StructuralMechanics3DKernel;

class StructuralMechanics3DControler: public StructuralMechanicsControler
{

public:
	StructuralMechanics3DControler(StructuralMechanicsLoadStepConfiguration &configuration);
	~StructuralMechanics3DControler();

	void dirichletIndices(std::vector<std::vector<esint> > &indices);
	void dirichletValues(std::vector<double> &values);

	void initData();
	void processSolution();

	void nextTime();
	void parametersChanged();

	void processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler);
	void processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler);

protected:
	StructuralMechanics3DKernel *_kernel;
};

}




#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS3D_CONTROLLER_H_ */
