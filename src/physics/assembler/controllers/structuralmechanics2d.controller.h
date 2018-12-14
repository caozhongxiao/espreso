
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS2D_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS2D_CONTROLLER_H_

#include "structuralmechanics.controller.h"

namespace espreso {

struct StructuralMechanics2DKernel;

class StructuralMechanics2DControler: public StructuralMechanicsControler
{

public:
	StructuralMechanics2DControler(StructuralMechanicsLoadStepConfiguration &configuration);
	~StructuralMechanics2DControler();

	void dirichletIndices(std::vector<std::vector<eslocal> > &indices);
	void dirichletValues(std::vector<double> &values);

	void analyticRegularization(size_t domain, bool ortogonalCluster);

	void initData();
	void processSolution();

	void nextTime();
	void parametersChanged();

	void processElements(Matrices matrices, InstanceFiller &filler);
	void processBoundary(Matrices matrices, size_t rindex, InstanceFiller &filler);

protected:
	StructuralMechanics2DKernel *_kernel;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS2D_CONTROLLER_H_ */
