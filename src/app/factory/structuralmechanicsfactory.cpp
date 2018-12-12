
#include "structuralmechanicsfactory.h"

#include "../../config/ecf/output.h"
#include "../../config/ecf/physics/structuralmechanics.h"

#include "../../physics/solver/timestep/linear.h"
#include "../../physics/solver/timestep/newtonraphson.h"
#include "../../physics/solver/loadstep/steadystate.h"

#include "../../physics/assembler/assembler.h"
#include "../../basis/logging/logging.h"
#include "../../physics/dataholder.h"

#include "../../physics/provider/distributedprovider.h"
#include "../../physics/provider/collectiveprovider.h"

using namespace espreso;

StructuralMechanicsFactory::StructuralMechanicsFactory(StructuralMechanicsConfiguration &configuration, ResultsSelectionConfiguration &propertiesConfiguration, Mesh *mesh)
: _configuration(configuration), _propertiesConfiguration(propertiesConfiguration), _bem(false)
{
	_instances.push_back(new DataHolder());

	switch (configuration.dimension) {
	case DIMENSION::D2:
		switch (configuration.load_steps_settings.at(1).solver) {
		case LoadStepConfiguration::SOLVER::FETI:
			_composer.push_back(new DomainsStructuralMechanics2D(
							*mesh, *_instances.front(), configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
			break;
		case LoadStepConfiguration::SOLVER::MULTIGRID:
			_composer.push_back(new GlobalStructuralMechanics3D(
							*mesh, *_instances.front(), configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
			break;
		}

		break;
	case DIMENSION::D3:
		switch (configuration.load_steps_settings.at(1).solver) {
		case LoadStepConfiguration::SOLVER::FETI:
			_composer.push_back(new DomainsStructuralMechanics3D(
							*mesh, *_instances.front(), configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
			break;
		case LoadStepConfiguration::SOLVER::MULTIGRID:
			_composer.push_back(new GlobalStructuralMechanics3D(
							*mesh, *_instances.front(), configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
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

LoadStepSolver* StructuralMechanicsFactory::getLoadStepSolver(size_t step, Mesh *mesh)
{
	StructuralMechanicsLoadStepConfiguration &settings = getLoadStepsSettings(step, _configuration.load_steps_settings);

	_linearSolvers.push_back(getLinearSolver(settings, _instances.front()));
	switch (_configuration.load_steps_settings.at(1).solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		_provider.push_back(new DistributedProvider(*_instances.front(), *_composer.front(), *mesh, *_linearSolvers.back()));
		break;
	case LoadStepConfiguration::SOLVER::MULTIGRID:
		_provider.push_back(new CollectiveProvider(*_instances.front(), *_composer.front(), *mesh, *_linearSolvers.back()));
		break;
	}

	switch (settings.mode) {
	case LoadStepConfiguration::MODE::LINEAR:
		_timeStepSolvers.push_back(new LinearTimeStep(*_provider.back()));
		break;
	case LoadStepConfiguration::MODE::NONLINEAR:
		if (_bem) {
			ESINFO(GLOBAL_ERROR) << "BEM discretization support only LINEAR STEADY STATE physics solver.";
		}
		switch (settings.nonlinear_solver.method) {
		case NonLinearSolverConfiguration::METHOD::NEWTON_RAPHSON:
		case NonLinearSolverConfiguration::METHOD::MODIFIED_NEWTON_RAPHSON:
			_timeStepSolvers.push_back(new NewtonRaphson(*_provider.back(), settings.nonlinear_solver));
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



