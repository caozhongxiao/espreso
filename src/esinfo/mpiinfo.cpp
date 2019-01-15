
#include "mpiinfo.h"

int espreso::info::mpi::MPIrank = 0;
int espreso::info::mpi::MPIsize = 1;
MPI_Comm espreso::info::mpi::MPICommunicator = MPI_COMM_WORLD;

void espreso::info::mpi::setMPI()
{
	int initialized;
	MPI_Initialized(&initialized);

	espreso::info::mpi::MPICommunicator = MPI_COMM_WORLD;
	if (initialized) {
		MPI_Comm_rank(espreso::info::mpi::MPICommunicator, &espreso::info::mpi::MPIrank);
		MPI_Comm_size(espreso::info::mpi::MPICommunicator, &espreso::info::mpi::MPIsize);
	}
}


