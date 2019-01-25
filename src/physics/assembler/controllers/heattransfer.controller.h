

#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER_CONTROLLER_H_

#include "controller.h"

namespace espreso {

struct HeatTransferLoadStepConfiguration;
struct NodeData;
struct ElementData;

class HeatTransferController: public Controller
{

public:
	void dirichletIndices(std::vector<std::vector<esint> > &indices);
	void dirichletValues(std::vector<double> &values);

	NodeData* solution();

protected:
	HeatTransferController(HeatTransferLoadStepConfiguration &configuration)
	: Controller(),
	  _configuration(configuration),
	  _temperature(NULL), _phaseChange(NULL), _latentHeat(NULL), _avgThickness(NULL),
	  _gradient(NULL), _flux(NULL), _motion(NULL), _thickness(NULL) {}

	HeatTransferLoadStepConfiguration &_configuration;

	struct BoundaryParameters {
		Parameter temperature, coordinate, thickness;
		Parameter htc, emissivity, externalTemperature, heatFlow, heatFlux;

		double regionArea;
	};

	Parameter _ntemperature, _ncoordinate, _nmotion, _nheat, _nthickness;
	std::vector<BoundaryParameters> _boundaries;

	NodeData *_temperature, *_phaseChange, *_latentHeat, *_avgThickness;
	ElementData *_gradient, *_flux, *_motion, *_thickness;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER_CONTROLLER_H_ */
