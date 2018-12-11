
#include "env.h"

int espreso::env::MPIrank = 0;
int espreso::env::MPIsize = 1;
MPI_Comm espreso::env::MPICommunicator = MPI_COMM_WORLD;

size_t espreso::env::MKL_NUM_THREADS = 1;
size_t espreso::env::OMP_NUM_THREADS = 1;
size_t espreso::env::SOLVER_NUM_THREADS = 1;
size_t espreso::env::PAR_NUM_THREADS = 1;

size_t espreso::env::verbose_level = 0;
size_t espreso::env::measure_level = 0;

using namespace espreso;

void env::setMPI()
{
	int initialized;
	MPI_Initialized(&initialized);

	env::MPICommunicator = MPI_COMM_WORLD;
	if (initialized) {
		MPI_Comm_rank(env::MPICommunicator, &env::MPIrank);
		MPI_Comm_size(env::MPICommunicator, &env::MPIsize);
	}
}




