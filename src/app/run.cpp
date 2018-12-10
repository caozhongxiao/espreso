
#include "mpi.h"

#include "factory/factory.h"

//#include "espreso.h"

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	espreso::ESPRESO::run(&argc, &argv);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

//	return espreso::run(&argc, &argv);
}


