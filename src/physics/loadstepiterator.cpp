
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

#include "assembler/provider/hypre/heattransfer.hypreprovider.h"
#include "assembler/provider/hypre/structuralmechanics.hypreprovider.h"
#include "assembler/provider/feti/heattransfer.fetiprovider.h"
#include "assembler/provider/feti/structuralmechanics2d.fetiprovider.h"
#include "assembler/provider/feti/structuralmechanics3d.fetiprovider.h"
#include "assembler/provider/direct/directprovider.h"

#include "basis/logging/logging.h"

#include "linearsolver/hypre/hypresolver.h"
#include "linearsolver/direct/directsolver.h"
#include "solver/generic/FETISolver.h"



using namespace espreso;

static Assembler* getAssembler(Assembler *previous, HeatTransferLoadStepConfiguration &loadStep, DIMENSION dimension)
{
	Assembler *current = NULL;
	switch (loadStep.solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		switch (dimension) {
		case DIMENSION::D2: current = new AssemblerInstance<HeatTransfer2DController, UniformNodesFETIComposer, HeatTransferFETIProvider, FETISolver>(previous, loadStep, loadStep.feti, 1); break;
		case DIMENSION::D3: current = new AssemblerInstance<HeatTransfer3DController, UniformNodesFETIComposer, HeatTransferFETIProvider, FETISolver>(previous, loadStep, loadStep.feti, 1); break;
		default: break;
		} break;
	case LoadStepConfiguration::SOLVER::HYPRE:
		switch (dimension) {
		case DIMENSION::D2: current = new AssemblerInstance<HeatTransfer2DController, UniformNodesComposer, HeatTransferHYPREProvider, HYPRESolver>(previous, loadStep, loadStep.hypre, 1); break;
		case DIMENSION::D3: current = new AssemblerInstance<HeatTransfer3DController, UniformNodesComposer, HeatTransferHYPREProvider, HYPRESolver>(previous, loadStep, loadStep.hypre, 1); break;
		default: break;
		} break;
	case LoadStepConfiguration::SOLVER::DIRECT:
		switch (dimension) {
		case DIMENSION::D2: current = new AssemblerInstance<HeatTransfer2DController, UniformNodesComposer, DirectProvider, DirectSolver>(previous, loadStep, 1); break;
		case DIMENSION::D3: current = new AssemblerInstance<HeatTransfer3DController, UniformNodesComposer, DirectProvider, DirectSolver>(previous, loadStep, 1); break;
		default: break;
		} break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented assembler.";
	}
	if (previous) {
		delete previous;
	}
	return current;
}

static Assembler* getAssembler(Assembler *previous, StructuralMechanicsLoadStepConfiguration &loadStep, DIMENSION dimension)
{
	Assembler *current = NULL;
	switch (loadStep.solver) {
	case LoadStepConfiguration::SOLVER::FETI:
		switch (dimension) {
		case DIMENSION::D2: current = new AssemblerInstance<StructuralMechanics2DController, UniformNodesFETIComposer, StructuralMechanics2DFETIProvider, FETISolver>(previous, loadStep, loadStep.feti, 2); break;
		case DIMENSION::D3: current = new AssemblerInstance<StructuralMechanics3DController, UniformNodesFETIComposer, StructuralMechanics3DFETIProvider, FETISolver>(previous, loadStep, loadStep.feti, 3); break;
		default: break;
		} break;
	case LoadStepConfiguration::SOLVER::HYPRE:
		switch (dimension) {
		case DIMENSION::D2: current = new AssemblerInstance<StructuralMechanics2DController, UniformNodesComposer, StructuralMechanicsHYPREProvider, HYPRESolver>(previous, loadStep, loadStep.hypre, 2); break;
		case DIMENSION::D3: current = new AssemblerInstance<StructuralMechanics3DController, UniformNodesComposer, StructuralMechanicsHYPREProvider, HYPRESolver>(previous, loadStep, loadStep.hypre, 3); break;
		default: break;
		} break;
	case LoadStepConfiguration::SOLVER::DIRECT:
		switch (dimension) {
		case DIMENSION::D2: current = new AssemblerInstance<StructuralMechanics2DController, UniformNodesComposer, DirectProvider, DirectSolver>(previous, loadStep, 1); break;
		case DIMENSION::D3: current = new AssemblerInstance<StructuralMechanics3DController, UniformNodesComposer, DirectProvider, DirectSolver>(previous, loadStep, 1); break;
		default: break;
		} break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented assembler.";
	}
	if (previous) {
		delete previous;
	}
	return current;
}

static TimeStepSolver* getTimeStepSolver(TimeStepSolver *previous, LoadStepConfiguration &loadStep, Assembler &assembler)
{
	TimeStepSolver* current = NULL;
	if (previous != NULL && !previous->hasSameMode(loadStep)) {
		delete previous;
		previous = NULL;
	}
	switch (loadStep.mode) {
	case LoadStepConfiguration::MODE::LINEAR: current = new LinearTimeStep(dynamic_cast<LinearTimeStep*>(previous), assembler); break;
	case LoadStepConfiguration::MODE::NONLINEAR: current = new NewtonRaphson(dynamic_cast<NewtonRaphson*>(previous), assembler, loadStep.nonlinear_solver); break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
	}
	if (previous) {
		delete previous;
	}
	return current;
}

static LoadStepSolver* getLoadStepSolver(LoadStepSolver *previous, HeatTransferLoadStepConfiguration &loadStep, Assembler &assembler, TimeStepSolver &timeStepSolver)
{
	LoadStepSolver* current = NULL;
	if (previous != NULL && !previous->hasSameType(loadStep)) {
		delete previous;
		previous = NULL;
	}
	switch (loadStep.type) {
	case LoadStepConfiguration::TYPE::STEADY_STATE:
		switch (loadStep.mode){
		case LoadStepConfiguration::MODE::LINEAR: current = new SteadyStateSolver(dynamic_cast<SteadyStateSolver*>(previous), assembler, timeStepSolver, loadStep.duration_time); break;
		case LoadStepConfiguration::MODE::NONLINEAR: current = new PseudoTimeStepping(dynamic_cast<PseudoTimeStepping*>(previous), assembler, timeStepSolver, loadStep.nonlinear_solver, loadStep.duration_time); break;
		} break;
	case LoadStepConfiguration::TYPE::TRANSIENT: current = new TransientFirstOrderImplicit(dynamic_cast<TransientFirstOrderImplicit*>(previous), assembler, timeStepSolver, loadStep.transient_solver, loadStep.duration_time); break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
	}
	if (previous) {
		delete previous;
	}
	return current;
}

static LoadStepSolver* getLoadStepSolver(LoadStepSolver *previous, StructuralMechanicsLoadStepConfiguration &loadStep, Assembler &assembler, TimeStepSolver &timeStepSolver)
{
	LoadStepSolver* current = NULL;
	if (previous != NULL && !previous->hasSameType(loadStep)) {
		delete previous;
		previous = NULL;
	}
	switch (loadStep.type) {
	case LoadStepConfiguration::TYPE::STEADY_STATE:
		switch (loadStep.mode){
		case LoadStepConfiguration::MODE::LINEAR: current = new SteadyStateSolver(dynamic_cast<SteadyStateSolver*>(previous), assembler, timeStepSolver, loadStep.duration_time); break;
		case LoadStepConfiguration::MODE::NONLINEAR: current = new PseudoTimeStepping(dynamic_cast<PseudoTimeStepping*>(previous), assembler, timeStepSolver, loadStep.nonlinear_solver, loadStep.duration_time); break;
		} break;
	case LoadStepConfiguration::TYPE::TRANSIENT: break;
//		return new TransientFirstOrderImplicit(assembler, timeStepSolver, loadStep.transient_solver, loadStep.duration_time);
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented requested SOLVER.";
	}
	if (previous) {
		delete previous;
	}
	return current;
}


LoadStepIterator::LoadStepIterator()
: _loadStepSolver(NULL), _timeStepSolver(NULL), _assembler(NULL)
{

}

LoadStepIterator::~LoadStepIterator()
{
	if (_loadStepSolver != NULL) { delete _loadStepSolver; }
	if (_timeStepSolver != NULL) { delete _timeStepSolver; }
	if (_assembler != NULL) { delete _assembler; }
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
	_assembler = getAssembler(_assembler, configuration.load_steps_settings.at(time::step + 1), configuration.dimension);
	_timeStepSolver = getTimeStepSolver(_timeStepSolver, configuration.load_steps_settings.at(time::step + 1), *_assembler);
	_loadStepSolver = getLoadStepSolver(_loadStepSolver, configuration.load_steps_settings.at(time::step + 1), *_assembler, *_timeStepSolver);

	_assembler->init();
	if (time::isInitial()) {
		info::mesh->storeMesh();
	}

	_loadStepSolver->run();

	return ++time::step < configuration.load_steps;
}


