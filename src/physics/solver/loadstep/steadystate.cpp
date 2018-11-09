
#include "../../solver/loadstep/steadystate.h"

#include "../../step.h"
#include "../../instance.h"
#include "../../solver/assembler.h"
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
	_assembler.updateStructuralMatrices(matrices);
	_assembler.updateGluingMatrices(matrices);
	return matrices;
}

void SteadyStateSolver::runNextTimeStep()
{
	_assembler.step.currentTime += _duration;
	_assembler.step.timeStep = _duration;
	processTimeStep();
}


void SteadyStateSolver::processTimeStep()
{
	_assembler.step.internalForceReduction = 1;
	_assembler.step.timeIntegrationConstantK = 1;
	_assembler.step.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(*this);

	_assembler.processSolution();
	_assembler.storeSolution();
}

