
#include "pseudotimestepping.h"
#include "physics/solver/timestep/timestepsolver.h"

#include "physics/dataholder.h"
#include "physics/assembler/assembler.h"

#include "globals/time.h"
#include "config/ecf/physics/physicssolver/nonlinearsolver.h"

using namespace espreso;

PseudoTimeStepping::PseudoTimeStepping(Assembler &assembler, TimeStepSolver &timeStepSolver, NonLinearSolverConfiguration &configuration, double duration)
: LoadStepSolver(assembler, timeStepSolver, duration), _configuration(configuration)
{

}

std::string PseudoTimeStepping::name()
{
	return "PSEUDO TIME STEPPING";
}

Matrices PseudoTimeStepping::updateStructuralMatrices(Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::f | Matrices::R | Matrices::Dirichlet);

	if (time::iteration) {
		updatedMatrices &= (Matrices::f | Matrices::Dirichlet | Matrices::R);
	}

	_assembler.assemble(updatedMatrices);
	_assembler.setDirichlet();
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
	_assembler.parameters.internalForceReduction = (double)(time::substep + 1) / _configuration.substeps;
	_assembler.parameters.timeIntegrationConstantK = 1;
	_assembler.parameters.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(*this);
}




