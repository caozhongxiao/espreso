
#include "lineartimesolver.h"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/assembler.h"
#include "physics/solver/loadstep/loadstepsolver.h"

#include "linearsolver/linearsolver.h"

#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

LinearTimeStep::LinearTimeStep(LinearTimeStep *previous, Assembler &assembler)
: TimeStepSolver(assembler)
{

}

bool LinearTimeStep::hasSameMode(const LoadStepConfiguration &configuration) const
{
	return configuration.mode == LoadStepConfiguration::MODE::LINEAR;
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
	_assembler.parametersChanged();
}




