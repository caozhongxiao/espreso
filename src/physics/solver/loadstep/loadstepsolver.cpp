
#include "../../solver/loadstep/loadstepsolver.h"

#include "../../../config/ecf/physics/heattransfer.h"
#include "../../../config/ecf/physics/structuralmechanics.h"

#include "../../../basis/logging/logging.h"
#include "../../../globals/time.h"
#include "../../provider/provider.h"
#include "../../solver/timestep/timestepsolver.h"

using namespace espreso;

LoadStepSolver::LoadStepSolver(const std::string &description, TimeStepSolver &timeStepSolver, double duration)
: _description(description), _timeStepSolver(timeStepSolver), _composer(timeStepSolver._composer), _duration(duration),
  _startTime(0), _precision(1e-8)
{

}

//LoadStepSolver* LoadStepSolver::create(const HeatTransferConfiguration &configuration)
//{
//	switch (configuration.dimension) {
//	case DIMENSION::D2:
//		break;
//	case DIMENSION::D3:
//		break;
//	}
//}
//
//LoadStepSolver* LoadStepSolver::create(const StructuralMechanicsConfiguration &configuration)
//{
//	switch (configuration.dimension) {
//	case DIMENSION::D2:
//		break;
//	case DIMENSION::D3:
//		break;
//	}
//}

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
	if (time::isInitial() == 0) {
		_composer.preprocessData();
	}
//	_composer.physics.setDirichlet();
	_composer.setRegularizationCallback();
	_composer.setB0Callback();
}

bool LoadStepSolver::hasNextTimeStep()
{
	return time::current + _precision < _startTime + _duration;
}

void LoadStepSolver::finalizeLoadStep()
{
	_composer.finalize();
}

void LoadStepSolver::run()
{
	ESINFO(PROGRESS1) << "Solve LOAD STEP " << time::step + 1 << ": " << description() << " with " << _timeStepSolver.description() << " time step(s).";

	_startTime = time::current;
	time::substep = 0;
	time::iteration = 0;

	initLoadStep();
	while (hasNextTimeStep()) {
		runNextTimeStep();
		ESINFO(PROGRESS1) << description() << " SOLVER: load step " << time::step + 1 << ", time step " << time::substep + 1 << " [" << time::current << "s] finished.";
		time::substep++;
		time::iteration = 0;
	}
	finalizeLoadStep();
}
