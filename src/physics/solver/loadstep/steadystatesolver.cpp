
#include "steadystatesolver.h"
#include "esinfo/timeinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "physics/assembler/dataholder.h"
#include "physics/solver/timestep/timestepsolver.h"
#include "physics/assembler/assembler.h"

#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

SteadyStateSolver::SteadyStateSolver(SteadyStateSolver *previous, Assembler &assembler, TimeStepSolver &timeStepSolver, double duration)
: LoadStepSolver(assembler, timeStepSolver, duration)
{

}

bool SteadyStateSolver::hasSameType(const LoadStepSolverConfiguration &configuration) const
{
	return
			configuration.type == LoadStepSolverConfiguration::TYPE::STEADY_STATE &&
			configuration.mode == LoadStepSolverConfiguration::MODE::LINEAR;
}

std::string SteadyStateSolver::name()
{
	return "STEADY STATE";
}

void SteadyStateSolver::setSolverParams()
{

}

Matrices SteadyStateSolver::updateStructuralMatrices(Matrices matrices)
{
	matrices &= (Matrices::K | Matrices::f | Matrices::R);
	_assembler.assemble(matrices);
	return matrices;
}

void SteadyStateSolver::runNextTimeStep()
{
	time::current += _duration;
	time::shift = _duration;
	_assembler.nextTime();

	processTimeStep();
}

void SteadyStateSolver::processTimeStep()
{
	_assembler.parameters.internalForceReduction = 1;
	_assembler.parameters.timeIntegrationConstantK = 1;
	_assembler.parameters.timeIntegrationConstantM = 0;

	eslog::solver("\n = ====================== STEADY STATE SOLVER ====================== =\n");
	eslog::solver(" =  LOAD STEP %2d, SUBSTEP %4d, TIME %10.6f, TIME STEP %8.6f  =\n", time::step + 1, time::substep + 1, time::current, time::shift);
	eslog::solver(" = ----------------------------------------------------------------- =\n");

	_timeStepSolver.solve(*this);
	_assembler.postProcess();

	eslog::solver(" = ================================================================= =\n");
	eslog::solver("                                             run time %12.3f s =\n\n", eslog::duration());
}

