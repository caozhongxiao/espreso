
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS3D_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS3D_CONTROLLER_H_

#include "structuralmechanics.controller.h"

namespace espreso {

struct StructuralMechanics3DKernel;

class StructuralMechanics3DController: public StructuralMechanicsController
{

public:
	StructuralMechanics3DController(StructuralMechanics3DController *previous, StructuralMechanicsLoadStepConfiguration &configuration);
	~StructuralMechanics3DController();

	const PhysicsConfiguration& configuration() const;

	void dirichletIndices(std::vector<std::vector<esint> > &indices);
	void dirichletValues(std::vector<double> &values);

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
