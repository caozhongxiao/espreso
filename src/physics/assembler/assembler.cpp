
#include "physics/assembler/dataholder.h"
#include "linearsolver/linearsolver.h"

#include "assembler.h"

using namespace espreso;

void Assembler::callsolve(Matrices matrices)
{
	_solver->solve(matrices);
}





