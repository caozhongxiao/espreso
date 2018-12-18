
#include "pseudotimestepping.h"
#include "../timestep/timestepsolver.h"

#include "../../dataholder.h"
#include "../../assembler/assembler.h"

#include "../../../globals/time.h"
#include "../../../config/ecf/physics/physicssolver/nonlinearsolver.h"

using namespace espreso;

PseudoTimeStepping::PseudoTimeStepping(Assembler &assembler, TimeStepSolver &timeStepSolver, NonLinearSolverConfiguration &configuration, double duration)
: LoadStepSolver(assembler, timeStepSolver, duration), _configuration(configuration)
{
//	_assembler.setRegularizationCallback();
//	_assembler.setB0Callback();
}

std::string PseudoTimeStepping::name()
{
	return "PSEUDO TRANSIENT";
}

Matrices PseudoTimeStepping::updateStructuralMatrices(Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::f | Matrices::R | Matrices::Dirichlet);

	if (time::iteration) {
		updatedMatrices &= (Matrices::f | Matrices::Dirichlet | Matrices::R);
	}

	_assembler.assemble(matrices);
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




