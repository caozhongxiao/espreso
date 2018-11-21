

#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLERS_HEATTRANSFER2D_CONTROLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLERS_HEATTRANSFER2D_CONTROLER_H_

#include "heattransfer.controler.h"

namespace espreso {

struct HeatTransfer2DKernel;

class HeatTransfer2DControler: public HeatTransferControler
{

public:
	HeatTransfer2DControler(
			Mesh &mesh, const Step &step,
			const HeatTransferGlobalSettings &gSettings,
			const HeatTransferStepSettings &sSettings,
			const HeatTransferOutputSettings &oSettings);
	~HeatTransfer2DControler();

	void initData();
	void updateData();

	void processElements(Matrices matrices, InstanceFiller &filler);
	void processBoundary(Matrices matrices, InstanceFiller &filler);

protected:
	void evaluate(const std::map<std::string, ECFExpression> &settings, tarray<double> &data);
	void evaluate(const std::map<std::string, ECFExpressionVector> &settings, tarray<double> &data);

	HeatTransfer2DKernel *_kernel;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLERS_HEATTRANSFER2D_CONTROLER_H_ */
