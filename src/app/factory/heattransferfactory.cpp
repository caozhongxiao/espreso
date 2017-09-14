
#include "../../config/ecf/physics/heattransfer.h"

#include "../../assembler/physics/laplacesteklovpoincare3d.h"

#include "../../assembler/physicssolver/assembler.h"
#include "../../assembler/physicssolver/timestep/linear.h"
#include "../../assembler/physicssolver/timestep/newtonraphson.h"
#include "../../assembler/physicssolver/loadstep/steadystate.h"
#include "../../assembler/physicssolver/loadstep/pseudotimestepping.h"
#include "../../assembler/physicssolver/loadstep/transientfirstorderimplicit.h"

#include "../../assembler/instance.h"
#include "../../assembler/physics/heattransfer2d.h"
#include "../../assembler/physics/heattransfer3d.h"
#include "../../mesh/structures/mesh.h"
#include "../../basis/logging/logging.h"
#include "heattransferfactory.h"

using namespace espreso;

HeatTransferFactory::HeatTransferFactory(const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration, Mesh *mesh)
: _configuration(configuration), _propertiesConfiguration(propertiesConfiguration), _bem(false)
{
	_instances.push_back(new Instance(*mesh));

	switch (configuration.discretization) {
	case DISCRETIZATION::FEM:
		switch (configuration.dimension) {
		case DIMENSION::D2:
			_physics.push_back(new HeatTransfer2D(mesh, _instances.front(), configuration, propertiesConfiguration));
			break;
		case DIMENSION::D3:
			_physics.push_back(new HeatTransfer3D(mesh, _instances.front(), configuration, propertiesConfiguration));
			break;
		}
		break;
	case DISCRETIZATION::BEM:
		_bem = true;
		switch (configuration.dimension) {
		case DIMENSION::D2:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: cannot solve HEAT TRANSFER 2D with BEM discretization.";
			break;
		case DIMENSION::D3:
			_physics.push_back(new LaplaceSteklovPoincare3D(mesh, _instances.front(), configuration, propertiesConfiguration));
			break;
		}
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown DISCRETIZATION for HeatTransfer3D";
	}
}

size_t HeatTransferFactory::loadSteps() const
{
	return _configuration.load_steps;
}

LoadStepSolver* HeatTransferFactory::getLoadStepSolver(size_t step, Mesh *mesh, Store *store)
{
	const HeatTransferLoadStepConfiguration &settings = getLoadStepsSettings(step, _configuration.load_steps_settings);

	_linearSolvers.push_back(getLinearSolver(settings, _instances.front()));
	_assemblers.push_back(new Assembler(*_instances.front(), *_physics.front(), *mesh, *store, *_linearSolvers.back()));

	switch (settings.mode) {
	case LoadStepConfiguration::MODE::LINEAR:
		_timeStepSolvers.push_back(new LinearTimeStep(*_assemblers.back()));
		break;
	case LoadStepConfiguration::MODE::NONLINEAR:
		if (_bem) {
			ESINFO(GLOBAL_ERROR) << "BEM discretization support only LINEAR STEADY STATE physics solver.";
		}
		switch (settings.nonlinear_solver.method) {
		case NonLinearSolverConfiguration::METHOD::NEWTON_RAPHSON:
		case NonLinearSolverConfiguration::METHOD::MODIFIED_NEWTON_RAPHSON:
			_timeStepSolvers.push_back(new NewtonRaphson(*_assemblers.back(), settings.nonlinear_solver));
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


