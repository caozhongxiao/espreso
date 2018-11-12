
#include "structuralmechanicsfactory.h"

#include "../../config/ecf/physics/structuralmechanics.h"
#include "../../physics/assembler/structuralmechanics2d.h"
#include "../../physics/assembler/structuralmechanics3d.h"
#include "../../physics/assembler/structuralmechanicstdnns3d.h"

#include "../../physics/composer/distributedcomposer.h"
#include "../../physics/solver/timestep/linear.h"
#include "../../physics/solver/timestep/newtonraphson.h"
#include "../../physics/solver/loadstep/steadystate.h"

#include "../../physics/instance.h"
#include "../../basis/logging/logging.h"

using namespace espreso;

StructuralMechanicsFactory::StructuralMechanicsFactory(Step *step, const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration, Mesh *mesh)
: _step(step), _configuration(configuration), _propertiesConfiguration(propertiesConfiguration), _bem(false)
{
	_instances.push_back(new Instance(*mesh));

	switch (configuration.dimension) {
	case DIMENSION::D2:
		_physics.push_back(new StructuralMechanics2D(mesh, _instances.front(), step, configuration, propertiesConfiguration));
		break;
	case DIMENSION::D3:
		switch (_configuration.assembler) {
		case ASSEMBLER::ELEMENTS:
			_physics.push_back(new StructuralMechanics3D(mesh, _instances.front(), step, configuration, propertiesConfiguration));
			break;
		case ASSEMBLER::FACES:
			_physics.push_back(new StructuralMechanicsTDNNS3D(mesh, _instances.front(), step, configuration, propertiesConfiguration));
			break;
		}

		break;
	default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid dimension.";
	}

	for (auto it = _configuration.discretization.begin(); it != _configuration.discretization.end(); ++it) {
		if (it->second == DISCRETIZATION::BEM) {
			_bem = true;
		}
	}
}

size_t StructuralMechanicsFactory::loadSteps() const
{
	return _configuration.load_steps;
}

LoadStepSolver* StructuralMechanicsFactory::getLoadStepSolver(size_t step, Mesh *mesh, ResultStore *store)
{
	const StructuralMechanicsLoadStepConfiguration &settings = getLoadStepsSettings(step, _configuration.load_steps_settings);

	_linearSolvers.push_back(getLinearSolver(settings, _instances.front()));
	_composers.push_back(new DistributedComposer(*_instances.front(), *_physics.front(), *mesh, *_step, *store, *_linearSolvers.back()));

	switch (settings.mode) {
	case LoadStepConfiguration::MODE::LINEAR:
		_timeStepSolvers.push_back(new LinearTimeStep(*_composers.back()));
		break;
	case LoadStepConfiguration::MODE::NONLINEAR:
		if (_bem) {
			ESINFO(GLOBAL_ERROR) << "BEM discretization support only LINEAR STEADY STATE physics solver.";
		}
		switch (settings.nonlinear_solver.method) {
		case NonLinearSolverConfiguration::METHOD::NEWTON_RAPHSON:
		case NonLinearSolverConfiguration::METHOD::MODIFIED_NEWTON_RAPHSON:
			_timeStepSolvers.push_back(new NewtonRaphson(*_composers.back(), settings.nonlinear_solver));
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented NONLINEAR SOLVER METHOD for LOAD STEP=" << step;
		}
		break;

	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented LOAD STEP solver MODE for LOAD STEP=" << step;
	}

	switch (settings.type) {
	case LoadStepConfiguration::TYPE::STEADY_STATE:
		_loadStepSolvers.push_back(new SteadyStateSolver(*_timeStepSolvers.back(), settings.duration_time));
		break;
	case LoadStepConfiguration::TYPE::TRANSIENT:
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented LOAD STEP solver TYPE for LOAD STEP=" << step;
	}

	return _loadStepSolvers.back();
}



