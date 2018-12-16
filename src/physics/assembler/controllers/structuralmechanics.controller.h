
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS_CONTROLLER_H_

#include "controller.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;
struct NodeData;
struct ElementData;

class StructuralMechanicsControler: public Controler
{

public:
	std::vector<double>& getSolutionStore();

protected:
	StructuralMechanicsControler(StructuralMechanicsLoadStepConfiguration &configuration)

	: Controler(),
	  _configuration(configuration),
	  _displacement(NULL), _avgThickness(NULL),
	  _thickness(NULL) {}

	const StructuralMechanicsLoadStepConfiguration &_configuration;

	struct BoundaryParameters {
		Parameter coordinate, thickness;
		Parameter normalPressure;
	};

	Parameter _ncoordinate, _ntemperature, _nInitialTemperature, _nthickness;
	std::vector<BoundaryParameters> _boundaries;

	NodeData *_displacement, *_avgThickness;
	ElementData *_thickness;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS_CONTROLLER_H_ */
