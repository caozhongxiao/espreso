
#include "../../solver/loadstep/loadstepsolver.h"

#include "../../step.h"
#include "../../../basis/logging/logging.h"
#include "../../assembler/physics.h"
#include "../../solver/assembler.h"
#include "../../solver/timestep/timestepsolver.h"

using namespace espreso;

LoadStepSolver::LoadStepSolver(const std::string &description, TimeStepSolver &timeStepSolver, double duration)
: _description(description), _timeStepSolver(timeStepSolver), _assembler(timeStepSolver._assembler), _duration(duration),
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
	if (_assembler.step.step == 0) {
		_assembler.preprocessData();
	}
	_assembler.physics.setDirichlet();
	_assembler.setRegularizationCallback();
	_assembler.setB0Callback();
}

bool LoadStepSolver::hasNextTimeStep()
{
	return _assembler.step.currentTime + _precision < _startTime + _duration;
}

void LoadStepSolver::finalizeLoadStep()
{
	_assembler.finalize();
}

void LoadStepSolver::run()
{
	ESINFO(PROGRESS1) << "Solve LOAD STEP " << _assembler.step.step + 1 << ": " << description() << " with " << _timeStepSolver.description() << " time step(s).";

	_startTime = _assembler.step.currentTime;
	_assembler.step.substep = 0;
	_assembler.step.iteration = 0;

	initLoadStep();
	while (hasNextTimeStep()) {
		runNextTimeStep();
		ESINFO(PROGRESS1) << description() << " SOLVER: load step " << _assembler.step.step + 1 << ", time step " << _assembler.step.substep + 1 << " [" << _assembler.step.currentTime << "s] finished.";
		_assembler.step.substep++;
		_assembler.step.iteration = 0;
	}
	finalizeLoadStep();
}
