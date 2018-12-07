
#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_STRUCTURALMECHANICS_CONTROLLER_H_

#include "controller.h"

namespace espreso {

struct StructuralMechanicsStepSettings;
struct StructuralMechanicsGlobalSettings;
struct StructuralMechanicsOutputSettings;
struct NodeData;
struct ElementData;

class StructuralMechanicsControler: public Controler
{

public:
	MatrixType getMatrixType(size_t domain) const;

	MatrixType getMatrixType() const;

	std::vector<double>& getSolutionStore();

protected:
	StructuralMechanicsControler(Mesh &mesh, const Step &step,
			const StructuralMechanicsGlobalSettings &gSettings,
			const StructuralMechanicsStepSettings &sSettings,
			const StructuralMechanicsOutputSettings &oSettings)

	: Controler(mesh, step),
	  _globalSettings(gSettings), _stepSettings(sSettings), _outputSettings(oSettings),
	  _displacement(NULL), _avgThickness(NULL),
	  _thickness(NULL) {}

	const StructuralMechanicsGlobalSettings &_globalSettings;
	const StructuralMechanicsStepSettings &_stepSettings;
	const StructuralMechanicsOutputSettings &_outputSettings;

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
