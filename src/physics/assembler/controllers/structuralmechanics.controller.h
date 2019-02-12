
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS_CONTROLLER_H_

#include "controller.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;
struct NodeData;
struct ElementData;

class StructuralMechanicsController: public Controller
{

public:
	NodeData* solution();

protected:
	StructuralMechanicsController(int dimension, StructuralMechanicsController *previous, StructuralMechanicsLoadStepConfiguration &configuration);

	StructuralMechanicsLoadStepConfiguration &_configuration;

	struct BoundaryParameters {
		Parameter coordinate, thickness;
		Parameter normalPressure;

		BoundaryParameters(int dimension)
		: coordinate(dimension), thickness(1),
		  normalPressure(1)
		{

		}
	};

	Parameter _kcoordinate, _ktemperature, _kinitialTemperature, _kacceleration, _kangularVelocity, _kthickness;
	std::vector<BoundaryParameters> _boundaries;

	NodeData *_ndisplacement;
	ElementData *_ethickness;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS_CONTROLLER_H_ */
