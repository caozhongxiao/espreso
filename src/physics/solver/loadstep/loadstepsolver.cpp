
#include "../../solver/loadstep/loadstepsolver.h"

#include "../../../basis/logging/logging.h"
#include "../../../globals/time.h"
#include "../../solver/timestep/timestepsolver.h"

using namespace espreso;

LoadStepSolver::LoadStepSolver(Assembler &assembler, TimeStepSolver &timeStepSolver, double duration)
: _assembler(assembler), _timeStepSolver(timeStepSolver), _duration(duration),
  _startTime(0), _precision(1e-8)
{

}

double LoadStepSolver::duration() const
{
	return _duration;
}

void LoadStepSolver::initLoadStep()
{
//	if (time::isInitial() == 0) {
//		_assembler.preprocessData();
//	}
////	_composer.physics.setDirichlet();
//	_assembler.setRegularizationCallback();
//	_assembler.setB0Callback();
}

bool LoadStepSolver::hasNextTimeStep()
{
	return time::current + _precision < _startTime + _duration;
}

void LoadStepSolver::finalizeLoadStep()
{
//	_assembler.finalize();
}

void LoadStepSolver::run()
{
	ESINFO(PROGRESS1) << "LOAD STEP: " << time::step << ", TYPE: " << name() << ", MODE: " << _timeStepSolver.name();

	_startTime = time::current;
	time::substep = 0;
	time::iteration = 0;

	initLoadStep();
	while (hasNextTimeStep()) {
		runNextTimeStep();
//		ESINFO(PROGRESS1) << description() << " SOLVER: load step " << time::step + 1 << ", time step " << time::substep + 1 << " [" << time::current << "s] finished.";
		time::substep++;
		time::iteration = 0;
	}
	finalizeLoadStep();
}
