
#include "mpiinfo.h"

int espreso::info::mpi::rank = 0;
int espreso::info::mpi::size = 1;
MPI_Comm espreso::info::mpi::comm = MPI_COMM_WORLD;

void espreso::info::mpi::set()
{
	int initialized;
	MPI_Initialized(&initialized);

	espreso::info::mpi::comm = MPI_COMM_WORLD;
	if (initialized) {
		MPI_Comm_rank(espreso::info::mpi::comm, &espreso::info::mpi::rank);
		MPI_Comm_size(espreso::info::mpi::comm, &espreso::info::mpi::size);
	}
}


