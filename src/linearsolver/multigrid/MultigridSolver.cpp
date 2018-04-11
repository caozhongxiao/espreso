/*
 * MultigridSolver.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: beh01
 */

#include "../../basis/utilities/utils.h"
#include "MultigridSolver.h"

#include "../../wrappers/hypre/hypre.h"

#include "../../config/ecf/environment.h"

using namespace espreso;

MultigridSolver::MultigridSolver(Instance *instance, const MultigridConfiguration &configuration)
: instance(instance),
  timeEvalMain("ESPRESO Multigrid Solver Overal Timing")
{
	this->configuration = configuration;
}

MultigridSolver::~MultigridSolver() {
	// TODO Auto-generated destructor stub
}


// make partial initialization according to updated matrices
void MultigridSolver::update(Matrices matrices)
{
	HYPRE::Solve(environment->MPIrank, environment->MPIsize, environment->MPICommunicator );
	ESINFO(GLOBAL_ERROR) << "Multigrid UPDATE not implemented.";
}


// run solver and store primal and dual solution
void MultigridSolver::solve()
{
	ESINFO(GLOBAL_ERROR) << "Multigrid SOLVE not implemented.";
}

void MultigridSolver::finalize()
{
	ESINFO(GLOBAL_ERROR) << "Multigrid FINALIZE not implemented.";
}
