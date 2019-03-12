
#include "loadstepsolver.h"
#include "esinfo/timeinfo.h"
#include "esinfo/eslog.h"

#include "physics/solver/timestep/timestepsolver.h"

using namespace espreso;

LoadStepSolver::LoadStepSolver(Assembler &assembler, TimeStepSolver &timeStepSolver, double duration)
: _assembler(assembler), _timeStepSolver(timeStepSolver), _duration(duration),
  _startTime(time::current), _precision(1e-8)
{

}

void LoadStepSolver::initLoadStep()
{

}

bool LoadStepSolver::hasNextTimeStep()
{
	return time::current + _precision < _startTime + _duration;
}

void LoadStepSolver::finalizeLoadStep()
{

}

void LoadStepSolver::run()
{
	eslog::start("SOLVER STARTED", "PHYSICS SOLVER");
	eslog::param("STEP", time::step + 1);
	eslog::param("TYPE", name().c_str());
	eslog::param("MODE", _timeStepSolver.name().c_str());
	eslog::ln();

	eslog::printsolverheader();

	_startTime = time::current;
	time::substep = 0;
	time::iteration = 0;

	initLoadStep();
	while (hasNextTimeStep()) {
		runNextTimeStep();
		time::substep++;
		time::iteration = 0;
	}
	finalizeLoadStep();

	eslog::endln("SOLVER FINISHED");
}
