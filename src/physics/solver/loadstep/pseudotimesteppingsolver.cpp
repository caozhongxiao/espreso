
#include "pseudotimesteppingsolver.h"
#include "esinfo/timeinfo.h"
#include "physics/assembler/dataholder.h"
#include "physics/solver/timestep/timestepsolver.h"

#include "physics/assembler/assembler.h"

#include "config/ecf/physics/physicssolver/nonlinearsolver.h"
#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

PseudoTimeStepping::PseudoTimeStepping(PseudoTimeStepping *previous, Assembler &assembler, TimeStepSolver &timeStepSolver, NonLinearSolverConfiguration &configuration, double duration)
: LoadStepSolver(assembler, timeStepSolver, duration), _configuration(configuration)
{

}

bool PseudoTimeStepping::hasSameType(const LoadStepConfiguration &configuration) const
{
	return
			configuration.type == LoadStepConfiguration::TYPE::STEADY_STATE &&
			configuration.mode == LoadStepConfiguration::MODE::NONLINEAR;
}

std::string PseudoTimeStepping::name()
{
	return "PSEUDO TIME STEPPING";
}

Matrices PseudoTimeStepping::updateStructuralMatrices(Matrices matrices)
{
	matrices &= (Matrices::K | Matrices::f | Matrices::R);
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
	_assembler.nextTime();

	processTimeStep();
}

void PseudoTimeStepping::processTimeStep()
{
	_assembler.parameters.internalForceReduction = (double)(time::substep + 1) / _configuration.substeps;
	_assembler.parameters.timeIntegrationConstantK = 1;
	_assembler.parameters.timeIntegrationConstantM = 0;

	_timeStepSolver.solve(*this);
	_assembler.postProcess();
}




