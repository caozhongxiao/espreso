

#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER2D_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER2D_CONTROLLER_H_

#include "heattransfer.controller.h"

namespace espreso {

struct HeatTransfer2DKernel;

class HeatTransfer2DControler: public HeatTransferControler
{

public:
	HeatTransfer2DControler(
			Mesh &mesh,
			const HeatTransferGlobalSettings &gSettings,
			const HeatTransferStepSettings &sSettings,
			const HeatTransferOutputSettings &oSettings);
	~HeatTransfer2DControler();

	void initData();
	void processSolution();

	void nextTime();
	void parametersChanged();

	void processElements(Matrices matrices, InstanceFiller &filler);
	void processBoundary(Matrices matrices, size_t rindex, InstanceFiller &filler);

protected:
	HeatTransfer2DKernel *_kernel;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER2D_CONTROLLER_H_ */
