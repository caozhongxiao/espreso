
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "physics/assembler/dataholder.h"
#include "linearsolver/linearsolver.h"

#include "assembler.h"

using namespace espreso;

void Assembler::callsolve(Matrices matrices)
{
	store("linsolver", matrices);
	_solver->solve(matrices);
	store("linsolver", Matrices::Solution);
}

void Assembler::store(const std::string &prefix, Matrices matrices)
{
	if (info::ecf->output.print_matrices) {
		eslog::debug("STORE MATRICES: %s\n", prefix.c_str());
		data()->store(prefix, matrices);
	}
}





