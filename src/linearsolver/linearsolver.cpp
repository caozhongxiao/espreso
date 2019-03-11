
#include "physics/assembler/dataholder.h"
#include "linearsolver.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"

#include "solver/generic/SparseMatrix.h"

#include <limits>
#include <cmath>

using namespace espreso;

LinearSolver::LinearSolver(DataHolder *data)
: _data(data)
{

}

LinearSolver::~LinearSolver()
{

}

void LinearSolver::solve(Matrices matrices)
{
	eslog::checkpointln("LINEAR SOLVER: DATA PREPARED");

	update(matrices);

	eslog::checkpointln("LINEAR SOLVER: DATA PREPROCESSED");

	solve();

	eslog::checkpointln("LINEAR SOLVER: SOLVED");

	double mmax = std::numeric_limits<double>::min(), gmax = std::numeric_limits<double>::min();
	#pragma omp parallel for reduction(max:mmax)
	for (size_t d = 0; d < _data->primalSolution.size(); d++) {
		for (size_t i = 0; i < _data->primalSolution[d].size(); i++) {
			mmax= std::max(mmax, std::fabs(_data->primalSolution[d][i]));
		}
	}

	MPI_Allreduce(&mmax, &gmax, 1, MPI_DOUBLE, MPI_MAX, info::mpi::comm);

	double dplaces = 1 / precision() / std::max(1., pow(10, std::ceil(std::log10(gmax))));

	#pragma omp parallel for
	for (size_t d = 0; d < _data->primalSolution.size(); d++) {
		for (size_t i = 0; i < _data->primalSolution[d].size(); i++) {
			_data->primalSolution[d][i] = std::trunc(dplaces * _data->primalSolution[d][i]) / dplaces;
		}
	}
}



