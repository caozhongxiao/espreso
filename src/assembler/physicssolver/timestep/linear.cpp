
#include "linear.h"

#include "../assembler.h"
#include "../loadstep/loadstepsolver.h"
#include "../../instance.h"

using namespace espreso;

LinearTimeStep::LinearTimeStep(Assembler &assembler)
: TimeStepSolver("LINEAR", assembler)
{

}

void LinearTimeStep::solve(LoadStepSolver &loadStepSolver)
{
	_assembler.solve(loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::B1));
}




