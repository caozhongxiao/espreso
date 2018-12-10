
#include "../../solver/loadstep/pseudotimestepping.h"

#include "../../../config/ecf/physics/physicssolver/nonlinearsolver.h"
#include "../../../globals/time.h"
#include "../../dataholder.h"
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
	double last = time::current;
	time::current += _duration / _configuration.substeps;
	if (time::current + _precision >= _startTime + _duration) {
		time::current = _startTime + _duration;
	}
	time::shift = time::current - last;
	processTimeStep();
}

void PseudoTimeStepping::processTimeStep()
{
//	_composer.step.internalForceReduction = (double)(time::substep + 1) / _configuration.substeps;
//	_composer.step.timeIntegrationConstantK = 1;
//	_composer.step.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(*this);
}




