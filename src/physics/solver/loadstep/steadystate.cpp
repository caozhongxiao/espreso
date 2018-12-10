
#include "../../solver/loadstep/steadystate.h"

#include "../../../globals/time.h"
#include "../../dataholder.h"
#include "../../provider/provider.h"
#include "../../solver/timestep/timestepsolver.h"

using namespace espreso;

SteadyStateSolver::SteadyStateSolver(TimeStepSolver &timeStepSolver, double duration)
: LoadStepSolver("STEADY STATE", timeStepSolver, duration)
{

}

Matrices SteadyStateSolver::updateStructuralMatrices(Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::f | Matrices::R | Matrices::B1 | Matrices::B1c | Matrices::B1duplicity);

	return reassembleStructuralMatrices(updatedMatrices);
}

Matrices SteadyStateSolver::reassembleStructuralMatrices(Matrices matrices)
{
	_composer.updateStructuralMatrices(matrices);
	_composer.updateGluingMatrices(matrices);
	return matrices;
}

void SteadyStateSolver::runNextTimeStep()
{
	time::current += _duration;
	time::shift = _duration;
	processTimeStep();
}


void SteadyStateSolver::processTimeStep()
{
//	_composer.step.internalForceReduction = 1;
//	_composer.step.timeIntegrationConstantK = 1;
//	_composer.step.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(*this);
}

