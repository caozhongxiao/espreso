
#include "linear.h"
#include "../loadstep/loadstepsolver.h"

#include "../../dataholder.h"
#include "../../assembler/assembler.h"

#include "../../../linearsolver/linearsolver.h"

using namespace espreso;

LinearTimeStep::LinearTimeStep(Assembler &assembler, LinearSolver &solver)
: TimeStepSolver(assembler, solver)
{

}

std::string LinearTimeStep::name()
{
	return "LINEAR";
}

void LinearTimeStep::solve(LoadStepSolver &loadStepSolver)
{
	_solver.solve(loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::Dirichlet));
	_assembler.postProcess();
}




