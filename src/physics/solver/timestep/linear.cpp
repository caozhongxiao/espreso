
#include "../../solver/timestep/linear.h"

#include "../../dataholder.h"
#include "../../provider/provider.h"
#include "../../solver/loadstep/loadstepsolver.h"

using namespace espreso;

LinearTimeStep::LinearTimeStep(Provider &composer)
: TimeStepSolver("LINEAR", composer)
{

}

void LinearTimeStep::solve(LoadStepSolver &loadStepSolver)
{
	_composer.nextTime();
	_composer.solve(loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::B1));
	_composer.parametersChanged();
	_composer.processSolution();
	_composer.storeSolution();
}




