

#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLERS_HEATTRANSFER_CONTROLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLERS_HEATTRANSFER_CONTROLER_H_

#include "controler.h"

namespace espreso {

struct HeatTransferStepSettings;
struct HeatTransferGlobalSettings;
struct HeatTransferOutputSettings;
struct NodeData;
struct ElementData;

class HeatTransferControler: public Controler
{

public:
	MatrixType getMatrixType() const;
	MatrixType getMatrixType(size_t domain) const;

protected:
	HeatTransferControler(Mesh &mesh, const Step &step,
			const HeatTransferGlobalSettings &gSettings,
			const HeatTransferStepSettings &sSettings,
			const HeatTransferOutputSettings &oSettings)

	: Controler(mesh, step),
	  _globalSettings(gSettings), _stepSettings(sSettings), _outputSettings(oSettings),
	  _temperature(NULL), _phaseChange(NULL), _latentHeat(NULL),
	  _gradient(NULL), _flux(NULL), _motion(NULL), _thickness(NULL) {}

	const HeatTransferGlobalSettings &_globalSettings;
	const HeatTransferStepSettings &_stepSettings;
	const HeatTransferOutputSettings &_outputSettings;

	Parameter _ntemperature, _ncoordinate, _nmotion, _nheat, _nthickness;

	NodeData *_temperature, *_phaseChange, *_latentHeat;
	ElementData *_gradient, *_flux, *_motion, *_thickness;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLERS_HEATTRANSFER_CONTROLER_H_ */
