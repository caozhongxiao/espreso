
#include "structuralmechanicsfactory.h"

#include "../../config/ecf/physics/structuralmechanics.h"

#include "../../assembler/physics/structuralmechanics2d.h"
#include "../../assembler/physics/structuralmechanics3d.h"
#include "../../assembler/physics/lamesteklovpoincare3d.h"

#include "../../assembler/physicssolver/assembler.h"
#include "../../assembler/physicssolver/timestep/linear.h"
#include "../../assembler/physicssolver/timestep/newtonraphson.h"
#include "../../assembler/physicssolver/loadstep/steadystate.h"

#include "../../assembler/instance.h"
#include "../../mesh/structures/mesh.h"
#include "../../basis/logging/logging.h"

using namespace espreso;

StructuralMechanicsFactory::StructuralMechanicsFactory(const StructuralMechanics2DConfiguration &configuration, Mesh *mesh)
: _configuration(configuration), _bem(false)
{
	_instances.push_back(new Instance(*mesh));
	_physics.push_back(new StructuralMechanics2D(mesh, _instances.front(), configuration));
}

StructuralMechanicsFactory::StructuralMechanicsFactory(const StructuralMechanics3DConfiguration &configuration, Mesh *mesh)
: _configuration(configuration), _bem(false)
{
	_instances.push_back(new Instance(*mesh));

	switch (configuration.discretization) {
	case DISCRETIZATION::FEM:
		_physics.push_back(new StructuralMechanics3D(mesh, _instances.front(), configuration));
		break;
	case DISCRETIZATION::BEM:
		_bem = true;
		_physics.push_back(new LameSteklovPoincare3D(mesh, _instances.front(), configuration));
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown DISCRETIZATION for StructuralMechanics3D";
	}
}

size_t StructuralMechanicsFactory::loadSteps() const
{
	return _configuration.physics_solver.load_steps;
}

LoadStepSolver* StructuralMechanicsFactory::getLoadStepSolver(size_t step, Mesh *mesh, Store *store)
{
	const StructuralMechanicsLoadStepsConfiguration &settings = getLoadStepsSettings(step, _configuration.physics_solver.load_steps_settings);

	_linearSolvers.push_back(getLinearSolver(settings, _instances.front()));
	_assemblers.push_back(new Assembler(*_instances.front(), *_physics.front(), *mesh, *store, *_linearSolvers.back()));

	switch (settings.mode) {
	case LoadStepsConfiguration::MODE::LINEAR:
		_timeStepSolvers.push_back(new LinearTimeStep(*_assemblers.back()));
		break;
	case LoadStepsConfiguration::MODE::NONLINEAR:
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
	case LoadStepsConfiguration::TYPE::STEADY_STATE:
		_loadStepSolvers.push_back(new SteadyStateSolver(*_timeStepSolvers.back(), settings.duration_time));
		break;
	case LoadStepsConfiguration::TYPE::TRANSIENT:
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented LOAD STEP solver TYPE for LOAD STEP=" << step;
	}

	return _loadStepSolvers.back();
}



