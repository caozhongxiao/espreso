
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER3D_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER3D_CONTROLLER_H_

#include "heattransfer.controller.h"

namespace espreso {

struct HeatTransfer3DKernel;
struct BEMData;

class HeatTransfer3DController: public HeatTransferController
{

public:
	HeatTransfer3DController(HeatTransferLoadStepConfiguration &configuration);
	~HeatTransfer3DController();

	const PhysicsConfiguration& configuration() const;

	void initData();
	void processSolution();

	void nextTime();
	void parametersChanged();

	void processBEMdomain(esint domain, double *values);
	void fillBEMInterior(esint domain, double *values);

	void processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler);
	void processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler);

protected:
	HeatTransfer3DKernel *_kernel;
	std::vector<BEMData*> _bem;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER3D_CONTROLLER_H_ */
