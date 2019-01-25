
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
	Matrices updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f);
	_assembler.setDirichlet(updatedMatrices);
	_assembler.solve(updatedMatrices | Matrices::Dirichlet);
}




