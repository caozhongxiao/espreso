
#include "loadstepiterator.h"
#include "dataholder.h"

#include "assembler/assembler.h"
#include "solver/timestep/linear.h"
#include "solver/timestep/newtonraphson.h"
#include "solver/loadstep/steadystate.h"
#include "solver/loadstep/pseudotimestepping.h"
#include "solver/loadstep/transientfirstorderimplicit.h"

#include "assembler/controllers/heattransfer2d.controller.h"
#include "assembler/controllers/heattransfer3d.controller.h"
#include "assembler/composer/global/uniformnodescomposer.h"
#include "assembler/provider/global/heattransfer.globalprovider.h"
#include "assembler/composer/feti/uniformnodesfeticomposer.h"
#include "assembler/provider/feti/heattransfer.fetiprovider.h"


#include "../globals/run.h"
#include "../globals/time.h"
#include "../basis/logging/logging.h"
#include "../config/ecf/root.h"

#include "../linearsolver/multigrid/multigrid.h"
#include "../solver/generic/FETISolver.h"

using namespace espreso;

static LinearSolver* getLinearSolver(LoadStepConfiguration &loadStep)
{
	switch (loadStep.solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		return new FETISolver(run::data, loadStep.feti);
	case LoadStepConfiguration::SOLVER::MULTIGRID:
		return new MultigridSolver(loadStep.multigrid);
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
		return NULL;
	}
}

static Assembler* getAssembler(HeatTransferLoadStepConfiguration &loadStep, DIMENSION dimension)
{
	switch (loadStep.solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		switch (dimension) {
		case DIMENSION::D2: return new AssemblerInstance<HeatTransfer2DControler, UniformNodesFETIComposer, HeatTransferFETIProvider>(loadStep, loadStep.feti);
		case DIMENSION::D3: return new AssemblerInstance<HeatTransfer3DControler, UniformNodesFETIComposer, HeatTransferFETIProvider>(loadStep, loadStep.feti);
		} break;
	case LoadStepConfiguration::SOLVER::MULTIGRID:
		switch (dimension) {
		case DIMENSION::D2: return new AssemblerInstance<HeatTransfer2DControler, UniformNodesComposer, HeatTransferGlobalProvider>(loadStep);
		case DIMENSION::D3: return new AssemblerInstance<HeatTransfer3DControler, UniformNodesComposer, HeatTransferGlobalProvider>(loadStep);
		}
	}
	ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
	return NULL;
}

static TimeStepSolver* getTimeStepSolver(LoadStepConfiguration &loadStep, Assembler &assembler, LinearSolver &solver)
{
	switch (loadStep.mode) {
	case LoadStepConfiguration::MODE::LINEAR:
		return new LinearTimeStep(assembler, solver);
	case LoadStepConfiguration::MODE::NONLINEAR:
		return new NewtonRaphson(assembler, solver, loadStep.nonlinear_solver);
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
		return NULL;
	}
}

static LoadStepSolver* getLoadStepSolver(HeatTransferLoadStepConfiguration &loadStep, Assembler &assembler, TimeStepSolver &timeStepSolver)
{
	switch (loadStep.type) {
	case LoadStepConfiguration::TYPE::STEADY_STATE:
		switch (loadStep.mode){
		case LoadStepConfiguration::MODE::LINEAR:
			return new SteadyStateSolver(assembler, timeStepSolver, loadStep.duration_time);
		case LoadStepConfiguration::MODE::NONLINEAR:
			return new PseudoTimeStepping(assembler, timeStepSolver, loadStep.nonlinear_solver, loadStep.duration_time);
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
		}
		return NULL;
	case LoadStepConfiguration::TYPE::TRANSIENT:
		return new TransientFirstOrderImplicit(assembler, timeStepSolver, loadStep.transient_solver, loadStep.duration_time);
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
		return NULL;
	}
}

LoadStepIterator::LoadStepIterator()
: _loadStepSolver(NULL), _timeStepSolver(NULL), _assembler(NULL), _linearSolver(NULL)
{
	run::data = new DataHolder();
}

bool LoadStepIterator::next()
{
	switch (run::ecf->physics) {
	case PHYSICS::HEAT_TRANSFER_2D:
		return next(run::ecf->heat_transfer_2d);
		break;
	case PHYSICS::HEAT_TRANSFER_3D:
		return next(run::ecf->heat_transfer_3d);
		break;
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
//		next();
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown physics.";
	}

	return false;
}

bool LoadStepIterator::next(HeatTransferConfiguration &configuration)
{
	if (time::step < configuration.load_steps) {
		_linearSolver = getLinearSolver(configuration.load_steps_settings.at(time::step + 1));
		_assembler = getAssembler(configuration.load_steps_settings.at(time::step + 1), configuration.dimension);
		_timeStepSolver = getTimeStepSolver(configuration.load_steps_settings.at(time::step + 1), *_assembler, *_linearSolver);
		_loadStepSolver = getLoadStepSolver(configuration.load_steps_settings.at(time::step + 1), *_assembler, *_timeStepSolver);

		_assembler->init();
		run::storeMesh();
		_loadStepSolver->run();
	}

	return ++time::step < configuration.load_steps;
}

bool LoadStepIterator::next(StructuralMechanicsConfiguration &configuration)
{
	if (time::step < configuration.load_steps) {



		++time::step;
	}

	return time::step < configuration.load_steps;
}


