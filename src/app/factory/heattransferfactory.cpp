
#include "heattransferfactory.h"

#include "../../config/ecf/output.h"
#include "../../config/ecf/physics/heattransfer.h"

#include "../../physics/solver/timestep/linear.h"
#include "../../physics/solver/timestep/newtonraphson.h"
#include "../../physics/solver/loadstep/steadystate.h"
#include "../../physics/solver/loadstep/pseudotimestepping.h"
#include "../../physics/solver/loadstep/transientfirstorderimplicit.h"

#include "../../physics/instance.h"
#include "../../physics/assembler/assembler.h"
#include "../../physics/assembler/heattransfer2d.h"
#include "../../physics/assembler/heattransfer3d.h"
#include "../../basis/logging/logging.h"

#include "../../physics/provider/distributedprovider.h"
#include "../../physics/provider/collectiveprovider.h"

using namespace espreso;

HeatTransferFactory::HeatTransferFactory(Step *step, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration, Mesh *mesh)
: _step(step), _configuration(configuration), _propertiesConfiguration(propertiesConfiguration), _bem(false)
{
	_instances.push_back(new Instance(*mesh));

	switch (configuration.dimension) {
	case DIMENSION::D2:
		_physics.push_back(new HeatTransfer2D(mesh, _instances.front(), step, configuration, propertiesConfiguration));
		switch (configuration.load_steps_settings.at(1).solver) {
		case LoadStepConfiguration::SOLVER::FETI:
			_composer.push_back(new DomainsHeatTransfer2D(
							*mesh, *_instances.front(), *step, configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
			break;
		case LoadStepConfiguration::SOLVER::HYPRE:
			_composer.push_back(new GlobalHeatTransfer2D(
							*mesh, *_instances.front(), *step, configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
			break;
		}

		break;
	case DIMENSION::D3:
		_physics.push_back(new HeatTransfer3D(mesh, _instances.front(), step, configuration, propertiesConfiguration));
		switch (configuration.load_steps_settings.at(1).solver) {
		case LoadStepConfiguration::SOLVER::FETI:
			_composer.push_back(new DomainsHeatTransfer3D(
							*mesh, *_instances.front(), *step, configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
			break;
		case LoadStepConfiguration::SOLVER::HYPRE:
			_composer.push_back(new GlobalHeatTransfer3D(
							*mesh, *_instances.front(), *step, configuration, configuration.load_steps_settings.at(1), propertiesConfiguration));
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

size_t HeatTransferFactory::loadSteps() const
{
	return _configuration.load_steps;
}

LoadStepSolver* HeatTransferFactory::getLoadStepSolver(size_t step, Mesh *mesh, ResultStore *store)
{
	const HeatTransferLoadStepConfiguration &settings = getLoadStepsSettings(step, _configuration.load_steps_settings);

	_linearSolvers.push_back(getLinearSolver(settings, _instances.front()));
	switch (_configuration.load_steps_settings.at(1).solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		_provider.push_back(new DistributedProvider(*_instances.front(), *_physics.front(), *_composer.front(), *mesh, *_step, *store, *_linearSolvers.back()));
		break;
	case LoadStepConfiguration::SOLVER::HYPRE:
		_provider.push_back(new CollectiveProvider(*_instances.front(), *_physics.front(), *_composer.front(), *mesh, *_step, *store, *_linearSolvers.back()));
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
		if (settings.mode == LoadStepConfiguration::MODE::NONLINEAR) {
			_loadStepSolvers.push_back(new PseudoTimeStepping(*_timeStepSolvers.back(), settings.nonlinear_solver, settings.duration_time));
		} else {
			_loadStepSolvers.push_back(new SteadyStateSolver(*_timeStepSolvers.back(), settings.duration_time));
		}
		break;
	case LoadStepConfiguration::TYPE::TRANSIENT:
		_loadStepSolvers.push_back(new TransientFirstOrderImplicit(*_timeStepSolvers.back(), settings.transient_solver, settings.duration_time));
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented LOAD STEP solver TYPE for LOAD STEP=" << step;
	}

	return _loadStepSolvers.back();
}


