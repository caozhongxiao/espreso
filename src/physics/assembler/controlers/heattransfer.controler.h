

#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLERS_HEATTRANSFER_CONTROLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLERS_HEATTRANSFER_CONTROLER_H_

#include "controler.h"

namespace espreso {

struct HeatTransferStepSettings;
struct HeatTransferGlobalSettings;

class HeatTransferControler: public Controler
{

public:
	MatrixType getMatrixType() const;
	MatrixType getMatrixType(size_t domain) const;

protected:
	HeatTransferControler(Mesh &mesh, const Step &step, const HeatTransferGlobalSettings &gSettings, const HeatTransferStepSettings &sSettings)
	: Controler(mesh, step), _globalSettings(gSettings), _stepSettings(sSettings) {}

	const HeatTransferGlobalSettings &_globalSettings;
	const HeatTransferStepSettings &_stepSettings;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLERS_HEATTRANSFER_CONTROLER_H_ */
