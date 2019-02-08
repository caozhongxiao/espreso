
#include "physics/assembler/dataholder.h"
#include "loadstepiterator.h"
#include "esinfo/time.h"
#include "esinfo/ecfinfo.h"

#include "assembler/assembler.h"
#include "solver/timestep/linear.h"
#include "solver/timestep/newtonraphson.h"
#include "solver/loadstep/steadystate.h"
#include "solver/loadstep/pseudotimestepping.h"
#include "solver/loadstep/transientfirstorderimplicit.h"

#include "assembler/controllers/heattransfer2d.controller.h"
#include "assembler/controllers/heattransfer3d.controller.h"
#include "assembler/controllers/structuralmechanics2d.controller.h"
#include "assembler/controllers/structuralmechanics3d.controller.h"

#include "assembler/composer/global/uniformnodescomposer.h"
#include "assembler/composer/feti/uniformnodesfeticomposer.h"

#include "assembler/provider/feti/heattransfer.fetiprovider.h"
#include "assembler/provider/feti/structuralmechanics2d.fetiprovider.h"
#include "assembler/provider/feti/structuralmechanics3d.fetiprovider.h"

#include "basis/logging/logging.h"

#include "linearsolver/hypre/hypresolver.h"
#include "solver/generic/FETISolver.h"
#include "assembler/provider/hypre/heattransfer.hypreprovider.h"
#include "assembler/provider/hypre/structuralmechanics.hypreprovider.h"


using namespace espreso;

static Assembler* getAssembler(HeatTransferLoadStepConfiguration &loadStep, DIMENSION dimension)
{
	switch (loadStep.solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		switch (dimension) {
		case DIMENSION::D2: return new AssemblerInstance<HeatTransfer2DController, UniformNodesFETIComposer, HeatTransferFETIProvider, FETISolver>(loadStep, loadStep.feti, 1);
		case DIMENSION::D3: return new AssemblerInstance<HeatTransfer3DController, UniformNodesFETIComposer, HeatTransferFETIProvider, FETISolver>(loadStep, loadStep.feti, 1);
		default: break;
		} break;
	case LoadStepConfiguration::SOLVER::HYPRE:
		switch (dimension) {
		case DIMENSION::D2: return new AssemblerInstance<HeatTransfer2DController, UniformNodesComposer, HeatTransferHYPREProvider, HYPRESolver>(loadStep, loadStep.hypre, 1);
		case DIMENSION::D3: return new AssemblerInstance<HeatTransfer3DController, UniformNodesComposer, HeatTransferHYPREProvider, HYPRESolver>(loadStep, loadStep.hypre, 1);
		default: break;
		}
	}
	ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
	return NULL;
}

static Assembler* getAssembler(StructuralMechanicsLoadStepConfiguration &loadStep, DIMENSION dimension)
{
	switch (loadStep.solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		switch (dimension) {
		case DIMENSION::D2: return new AssemblerInstance<StructuralMechanics2DController, UniformNodesFETIComposer, StructuralMechanics2DFETIProvider, FETISolver>(loadStep, loadStep.feti, 2);
		case DIMENSION::D3: return new AssemblerInstance<StructuralMechanics3DController, UniformNodesFETIComposer, StructuralMechanics3DFETIProvider, FETISolver>(loadStep, loadStep.feti, 3);
		default: break;
		} break;
	case LoadStepConfiguration::SOLVER::HYPRE:
		switch (dimension) {
		case DIMENSION::D2: return new AssemblerInstance<StructuralMechanics2DController, UniformNodesComposer, StructuralMechanicsHYPREProvider, HYPRESolver>(loadStep, loadStep.hypre, 2);
		case DIMENSION::D3: return new AssemblerInstance<StructuralMechanics3DController, UniformNodesComposer, StructuralMechanicsHYPREProvider, HYPRESolver>(loadStep, loadStep.hypre, 3);
		default: break;
		}
	}
	ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
	return NULL;
}

static TimeStepSolver* getTimeStepSolver(LoadStepConfiguration &loadStep, Assembler &assembler, LinearSolver &solver)
{
	switch (loadStep.mode) {
	case LoadStepConfiguration::MODE::LINEAR:
		return new LinearTimeStep(assembler);
	case LoadStepConfiguration::MODE::NONLINEAR:
		return new NewtonRaphson(assembler, loadStep.nonlinear_solver);
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
		case LoadStepConfiguration::MODE::LINEAR: return new SteadyStateSolver(assembler, timeStepSolver, loadStep.duration_time);
		case LoadStepConfiguration::MODE::NONLINEAR: return new PseudoTimeStepping(assembler, timeStepSolver, loadStep.nonlinear_solver, loadStep.duration_time);
		} break;
	case LoadStepConfiguration::TYPE::TRANSIENT:
		return new TransientFirstOrderImplicit(assembler, timeStepSolver, loadStep.transient_solver, loadStep.duration_time);
	}
	ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
	return NULL;
}

static LoadStepSolver* getLoadStepSolver(StructuralMechanicsLoadStepConfiguration &loadStep, Assembler &assembler, TimeStepSolver &timeStepSolver)
{
	switch (loadStep.type) {
	case LoadStepConfiguration::TYPE::STEADY_STATE:
		switch (loadStep.mode){
		case LoadStepConfiguration::MODE::LINEAR: return new SteadyStateSolver(assembler, timeStepSolver, loadStep.duration_time);
		case LoadStepConfiguration::MODE::NONLINEAR: return new PseudoTimeStepping(assembler, timeStepSolver, loadStep.nonlinear_solver, loadStep.duration_time);
		} break;
	case LoadStepConfiguration::TYPE::TRANSIENT: break;
//		return new TransientFirstOrderImplicit(assembler, timeStepSolver, loadStep.transient_solver, loadStep.duration_time);
	}
	ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
	return NULL;
}


LoadStepIterator::LoadStepIterator()
: _loadStepSolver(NULL), _timeStepSolver(NULL), _assembler(NULL), _linearSolver(NULL)
{

}

LoadStepIterator::~LoadStepIterator()
{
	if (_loadStepSolver != NULL) { delete _loadStepSolver; }
	if (_timeStepSolver != NULL) { delete _timeStepSolver; }
	if (_assembler != NULL) { delete _assembler; }
	if (_linearSolver != NULL) { delete _linearSolver; }
}

bool LoadStepIterator::next()
{
	switch (info::ecf->physics) {
	case PHYSICS::HEAT_TRANSFER_2D:
		return next(info::ecf->heat_transfer_2d);
	case PHYSICS::HEAT_TRANSFER_3D:
		return next(info::ecf->heat_transfer_3d);
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
		return next(info::ecf->structural_mechanics_2d);
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
		return next(info::ecf->structural_mechanics_3d);
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown physics.";
	}

	return false;
}

template <typename TPhysics>
bool LoadStepIterator::next(TPhysics &configuration)
{
	if (time::isInitial()) {
		_assembler = getAssembler(configuration.load_steps_settings.at(time::step + 1), configuration.dimension);
		_timeStepSolver = getTimeStepSolver(configuration.load_steps_settings.at(time::step + 1), *_assembler, *_linearSolver);
		_loadStepSolver = getLoadStepSolver(configuration.load_steps_settings.at(time::step + 1), *_assembler, *_timeStepSolver);

		_assembler->init();
		info::mesh->storeMesh();
		_loadStepSolver->run();
	} else if (time::step < configuration.load_steps) {
//		_assembler = getAssembler(configuration.load_steps_settings.at(time::step + 1), configuration.dimension);
//		_timeStepSolver = getTimeStepSolver(configuration.load_steps_settings.at(time::step + 1), *_assembler, *_linearSolver);
//		_loadStepSolver = getLoadStepSolver(configuration.load_steps_settings.at(time::step + 1), *_assembler, *_timeStepSolver);
//
//		_loadStepSolver->run();
	}

	return ++time::step < configuration.load_steps;
}


