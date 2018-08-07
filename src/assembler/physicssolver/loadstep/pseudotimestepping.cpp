
#include "pseudotimestepping.h"

#include "../assembler.h"
#include "../timestep/timestepsolver.h"
#include "../../step.h"
#include "../../instance.h"
#include "../../../config/ecf/physics/physicssolver/nonlinearsolver.h"


using namespace espreso;

PseudoTimeStepping::PseudoTimeStepping(TimeStepSolver &timeStepSolver, const NonLinearSolverConfiguration &configuration, double duration)
: LoadStepSolver("PSEUDO TIME STEPS", timeStepSolver, duration), _configuration(configuration)
{
	_assembler.setRegularizationCallback();
	_assembler.setB0Callback();
}

Matrices PseudoTimeStepping::updateStructuralMatrices(Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::f | Matrices::R | Matrices::B1 | Matrices::B1c | Matrices::B1duplicity);

//	if (_assembler.step.iteration) {
//		updatedMatrices &= (Matrices::f | Matrices::B1c | Matrices::R);
//	}

	return reassembleStructuralMatrices(updatedMatrices);
}

Matrices PseudoTimeStepping::reassembleStructuralMatrices(Matrices matrices)
{
	_assembler.updateStructuralMatrices(matrices);
	_assembler.updateGluingMatrices(matrices);
	return matrices;
}

void PseudoTimeStepping::runNextTimeStep()
{
	double last = _assembler.step.currentTime;
	_assembler.step.currentTime += _duration / _configuration.substeps;
	if (_assembler.step.currentTime + _precision >= _startTime + _duration) {
		_assembler.step.currentTime = _startTime + _duration;
	}
	_assembler.step.timeStep = _assembler.step.currentTime - last;
	processTimeStep();
}

void PseudoTimeStepping::processTimeStep()
{
	_assembler.step.internalForceReduction = (double)(_assembler.step.substep + 1) / _configuration.substeps;
	_assembler.step.timeIntegrationConstantK = 1;
	_assembler.step.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(*this);

	_assembler.processSolution();
	_assembler.storeSolution();
}




