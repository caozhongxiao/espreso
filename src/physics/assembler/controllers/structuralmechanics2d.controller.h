
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS2D_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS2D_CONTROLLER_H_

#include "structuralmechanics.controller.h"

namespace espreso {

struct StructuralMechanics2DKernel;

class StructuralMechanics2DController: public StructuralMechanicsController
{

public:
	StructuralMechanics2DController(StructuralMechanicsLoadStepConfiguration &configuration);
	~StructuralMechanics2DController();

	const PhysicsConfiguration& configuration() const;

	void dirichletIndices(std::vector<std::vector<esint> > &indices);
	void dirichletValues(std::vector<double> &values);

	void initData();
	void processSolution();

	void nextTime();
	void parametersChanged();

	void processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler);
	void processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler);

protected:
	StructuralMechanics2DKernel *_kernel;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS2D_CONTROLLER_H_ */
