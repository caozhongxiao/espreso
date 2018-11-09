
#include "../../solver/timestep/linear.h"

#include "../../instance.h"
#include "../../solver/assembler.h"
#include "../../solver/loadstep/loadstepsolver.h"

using namespace espreso;

LinearTimeStep::LinearTimeStep(Assembler &assembler)
: TimeStepSolver("LINEAR", assembler)
{

}

void LinearTimeStep::solve(LoadStepSolver &loadStepSolver)
{
	_assembler.solve(loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::B1));
}




