
#include "linear.h"
#include "physics/assembler/dataholder.h"
#include "physics/assembler/assembler.h"
#include "physics/solver/loadstep/loadstepsolver.h"

#include "linearsolver/linearsolver.h"

using namespace espreso;

LinearTimeStep::LinearTimeStep(Assembler &assembler)
: TimeStepSolver(assembler)
{

}

std::string LinearTimeStep::name()
{
	return "LINEAR";
}

void LinearTimeStep::solve(LoadStepSolver &loadStepSolver)
{
	_assembler.solve(loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::Dirichlet));
}




