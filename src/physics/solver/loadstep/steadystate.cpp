
#include "physics/assembler/dataholder.h"
#include "esinfo/time.h"
#include "esinfo/meshinfo.h"
#include "steadystate.h"
#include "physics/solver/timestep/timestepsolver.h"
#include "physics/assembler/assembler.h"


using namespace espreso;

SteadyStateSolver::SteadyStateSolver(Assembler &assembler, TimeStepSolver &timeStepSolver, double duration)
: LoadStepSolver(assembler, timeStepSolver, duration)
{

}

std::string SteadyStateSolver::name()
{
	return "STEADY STATE";
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

	_timeStepSolver.solve(*this);
	_assembler.postProcess();
}

