
#include "../../solver/loadstep/steadystate.h"
#include "../../solver/timestep/timestepsolver.h"
#include "../../dataholder.h"
#include "../../assembler/assembler.h"

#include "../../../globals/time.h"


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
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::f | Matrices::R | Matrices::B1 | Matrices::B1c | Matrices::B1duplicity);

	return reassembleStructuralMatrices(updatedMatrices);
}

Matrices SteadyStateSolver::reassembleStructuralMatrices(Matrices matrices)
{
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
//	_composer.step.internalForceReduction = 1;
//	_composer.step.timeIntegrationConstantK = 1;
//	_composer.step.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(*this);
}

