
#include "../../solver/loadstep/pseudotimestepping.h"

#include "../../step.h"
#include "../../instance.h"
#include "../../../config/ecf/physics/physicssolver/nonlinearsolver.h"
#include "../../provider/provider.h"
#include "../../solver/timestep/timestepsolver.h"


using namespace espreso;

PseudoTimeStepping::PseudoTimeStepping(TimeStepSolver &timeStepSolver, const NonLinearSolverConfiguration &configuration, double duration)
: LoadStepSolver("PSEUDO TIME STEPS", timeStepSolver, duration), _configuration(configuration)
{
	_composer.setRegularizationCallback();
	_composer.setB0Callback();
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
	_composer.updateStructuralMatrices(matrices);
	_composer.updateGluingMatrices(matrices);
	return matrices;
}

void PseudoTimeStepping::runNextTimeStep()
{
	double last = _composer.step.currentTime;
	_composer.step.currentTime += _duration / _configuration.substeps;
	if (_composer.step.currentTime + _precision >= _startTime + _duration) {
		_composer.step.currentTime = _startTime + _duration;
	}
	_composer.step.timeStep = _composer.step.currentTime - last;
	processTimeStep();
}

void PseudoTimeStepping::processTimeStep()
{
	_composer.step.internalForceReduction = (double)(_composer.step.substep + 1) / _configuration.substeps;
	_composer.step.timeIntegrationConstantK = 1;
	_composer.step.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(*this);
}




