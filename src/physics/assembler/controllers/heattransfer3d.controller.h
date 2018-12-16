
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER3D_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER3D_CONTROLLER_H_

#include "heattransfer.controller.h"

namespace espreso {

struct HeatTransfer3DKernel;

class HeatTransfer3DControler: public HeatTransferControler
{

public:
	HeatTransfer3DControler(HeatTransferLoadStepConfiguration &configuration);
	~HeatTransfer3DControler();

	void initData();
	void processSolution();

	void nextTime();
	void parametersChanged();

	void processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler);
	void processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler);

protected:
	HeatTransfer3DKernel *_kernel;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER3D_CONTROLLER_H_ */
