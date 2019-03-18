
#include "lineartimesolver.h"

#include "esinfo/eslog.hpp"
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

bool LinearTimeStep::hasSameMode(const LoadStepSolverConfiguration &configuration) const
{
	return configuration.mode == LoadStepSolverConfiguration::MODE::LINEAR;
}

std::string LinearTimeStep::name()
{
	return "LINEAR";
}

void LinearTimeStep::setSolverParams()
{

}

void LinearTimeStep::solve(LoadStepSolver &loadStepSolver)
{
	Matrices updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f);
	_assembler.setDirichlet(updatedMatrices);

	eslog::solver("   - LINEAR TIME STEP        REASSEMBLED MATRICES :: %c, %c, %c, %c, %c -\n",
			(updatedMatrices & Matrices::K) ? 'K' : ' ',
			(updatedMatrices & Matrices::M) ? 'M' : ' ',
			(updatedMatrices & Matrices::C) ? 'C' : ' ',
			(updatedMatrices & Matrices::R) ? 'R' : ' ',
			(updatedMatrices & Matrices::f) ? 'f' : ' ');

	_assembler.solve(updatedMatrices | Matrices::Dirichlet);
	_assembler.parametersChanged();
}




