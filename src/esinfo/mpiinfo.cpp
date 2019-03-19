
#include "mpiinfo.h"
#include "eslog.hpp"

int espreso::info::mpi::rank = 0;
int espreso::info::mpi::size = 1;
MPI_Comm espreso::info::mpi::comm = MPI_COMM_WORLD;

int espreso::info::mpi::irank = 0;
int espreso::info::mpi::isize = 1;
MPI_Comm espreso::info::mpi::icomm = MPI_COMM_SELF;

int espreso::info::mpi::grank = 0;
int espreso::info::mpi::gsize = 1;
MPI_Comm espreso::info::mpi::gcomm = MPI_COMM_WORLD;

using namespace espreso::info;


void mpi::set()
{
	int initialized;
	MPI_Initialized(&initialized);

	mpi::comm = MPI_COMM_WORLD;
	if (initialized) {
		MPI_Comm_rank(mpi::comm, &mpi::rank);
		MPI_Comm_size(mpi::comm, &mpi::size);
	}

	mpi::irank = mpi::rank;
	mpi::isize = mpi::size;
	mpi::icomm = mpi::comm;

	mpi::grank = mpi::rank;
	mpi::gsize = mpi::size;
	mpi::gcomm = mpi::comm;
}

bool mpi::divide(int meshDuplication)
{
	if (meshDuplication == 1) {
		return true;
	}

	if (espreso::info::mpi::size % meshDuplication != 0) {
		return false;
	}

	int color = mpi::rank / (mpi::size / meshDuplication);

	MPI_Comm_split(mpi::gcomm, color, mpi::grank, &mpi::comm);
	MPI_Comm_rank(mpi::comm, &mpi::rank);
	MPI_Comm_size(mpi::comm, &mpi::size);

	MPI_Comm_split(mpi::gcomm, mpi::rank, mpi::grank, &mpi::icomm);
	MPI_Comm_rank(mpi::icomm, &mpi::irank);
	MPI_Comm_size(mpi::icomm, &mpi::isize);

	return true;
}

void mpi::finish()
{
	if (mpi::isize > 1) {
		MPI_Comm_free(&mpi::comm);
		MPI_Comm_free(&mpi::icomm);
	}
}


