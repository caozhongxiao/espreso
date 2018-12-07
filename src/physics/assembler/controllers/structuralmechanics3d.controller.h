
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS3D_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS3D_CONTROLLER_H_

#include "structuralmechanics.controller.h"

namespace espreso {

struct StructuralMechanics3DKernel;

class StructuralMechanics3DControler: public StructuralMechanicsControler
{

public:
	StructuralMechanics3DControler(
			Mesh &mesh, const Step &step,
			const StructuralMechanicsGlobalSettings &gSettings,
			const StructuralMechanicsStepSettings &sSettings,
			const StructuralMechanicsOutputSettings &oSettings);
	~StructuralMechanics3DControler();

	void analyticRegularization(size_t domain, bool ortogonalCluster);

	void dirichletIndices(std::vector<std::vector<eslocal> > &indices);
	void dirichletValues(std::vector<double> &values);

	void initData();
	void processSolution();

	void nextTime();
	void parametersChanged();

	void processElements(Matrices matrices, InstanceFiller &filler);
	void processBoundary(Matrices matrices, size_t rindex, InstanceFiller &filler);

protected:
	StructuralMechanics3DKernel *_kernel;
};

}




#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS3D_CONTROLLER_H_ */
