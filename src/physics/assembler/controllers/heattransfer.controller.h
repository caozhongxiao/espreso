

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
	HeatTransferController(int dimension, HeatTransferController *previous, HeatTransferLoadStepConfiguration &configuration);

	HeatTransferLoadStepConfiguration &_configuration;

	struct BoundaryParameters {
		Parameter temperature, coordinate, thickness;
		Parameter htc, emissivity, externalTemperature, heatFlow, heatFlux;

		double regionArea;

		BoundaryParameters(int dimension)
		: temperature(1), coordinate(dimension), thickness(1),
		  htc(1), emissivity(1), externalTemperature(1), heatFlow(1), heatFlux(1),
		  regionArea(0)
		{

		}
	};

	Parameter _ktemperature, _kcoordinate, _kmotion, _kheat, _kthickness;
	std::vector<BoundaryParameters> _boundaries;

	NodeData *_ntemperature;
	ElementData *_egradient, *_eflux, *_emotion, *_ethickness, *_ephaseChange, *_elatentHeat;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER_CONTROLLER_H_ */
