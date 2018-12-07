
#include "../../solver/loadstep/loadstepsolver.h"

#include "../../step.h"
#include "../../../basis/logging/logging.h"
#include "../../provider/provider.h"
#include "../../solver/timestep/timestepsolver.h"

using namespace espreso;

LoadStepSolver::LoadStepSolver(const std::string &description, TimeStepSolver &timeStepSolver, double duration)
: _description(description), _timeStepSolver(timeStepSolver), _composer(timeStepSolver._composer), _duration(duration),
  _startTime(0), _precision(1e-8)
{

}

std::string LoadStepSolver::description() const
{
	return _description;
}

double LoadStepSolver::duration() const
{
	return _duration;
}

void LoadStepSolver::initLoadStep()
{
	if (_composer.step.step == 0) {
		_composer.preprocessData();
	}
//	_composer.physics.setDirichlet();
	_composer.setRegularizationCallback();
	_composer.setB0Callback();
}

bool LoadStepSolver::hasNextTimeStep()
{
	return _composer.step.currentTime + _precision < _startTime + _duration;
}

void LoadStepSolver::finalizeLoadStep()
{
	_composer.finalize();
}

void LoadStepSolver::run()
{
	ESINFO(PROGRESS1) << "Solve LOAD STEP " << _composer.step.step + 1 << ": " << description() << " with " << _timeStepSolver.description() << " time step(s).";

	_startTime = _composer.step.currentTime;
	_composer.step.substep = 0;
	_composer.step.iteration = 0;

	initLoadStep();
	while (hasNextTimeStep()) {
		runNextTimeStep();
		ESINFO(PROGRESS1) << description() << " SOLVER: load step " << _composer.step.step + 1 << ", time step " << _composer.step.substep + 1 << " [" << _composer.step.currentTime << "s] finished.";
		_composer.step.substep++;
		_composer.step.iteration = 0;
	}
	finalizeLoadStep();
}
