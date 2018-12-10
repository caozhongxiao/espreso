

#ifndef SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER_CONTROLLER_H_
#define SRC_PHYSICS_ASSEMBLER_CONTROLLERS_HEATTRANSFER_CONTROLLER_H_

#include "controller.h"

namespace espreso {

struct HeatTransferStepSettings;
struct HeatTransferGlobalSettings;
struct HeatTransferOutputSettings;
struct NodeData;
struct ElementData;

class HeatTransferControler: public Controler
{

public:
	MatrixType getMatrixType(size_t domain) const;
	void analyticRegularization(size_t domain, bool ortogonalCluster);

	MatrixType getMatrixType() const;

	void dirichletIndices(std::vector<std::vector<eslocal> > &indices);
	void dirichletValues(std::vector<double> &values);

	std::vector<double>& getSolutionStore();

protected:
	HeatTransferControler(Mesh &mesh,
			const HeatTransferGlobalSettings &gSettings,
			const HeatTransferStepSettings &sSettings,
			const HeatTransferOutputSettings &oSettings)

	: Controler(mesh),
	  _globalSettings(gSettings), _stepSettings(sSettings), _outputSettings(oSettings),
	  _temperature(NULL), _phaseChange(NULL), _latentHeat(NULL), _avgThickness(NULL),
	  _gradient(NULL), _flux(NULL), _motion(NULL), _thickness(NULL) {}

	const HeatTransferGlobalSettings &_globalSettings;
	const HeatTransferStepSettings &_stepSettings;
	const HeatTransferOutputSettings &_outputSettings;

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
