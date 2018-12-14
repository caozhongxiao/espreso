
#include "../../solver/timestep/linear.h"
#include "../../solver/loadstep/loadstepsolver.h"

#include "../../dataholder.h"
#include "../../assembler/assembler.h"

#include "../../../linearsolver/linearsolver.h"

#include "../../../globals/run.h"
#include "../../../mesh/mesh.h"
#include "../../../output/result/resultstore.h"

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
	_solver.update(loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::B1));
	_solver.solve();
	_assembler.postProcess();
	run::mesh->store->updateSolution();
}




