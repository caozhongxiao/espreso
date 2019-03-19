
#ifndef SRC_ESINFO_MPIINFO_H_
#define SRC_ESINFO_MPIINFO_H_

#include "mpi.h"

namespace espreso {
namespace info {
namespace mpi {
	// instance MPI communication settings
	extern int rank;
	extern int size;
	extern MPI_Comm comm;

	// inter-instances MPI communication settings
	extern int irank;
	extern int isize;
	extern MPI_Comm icomm;

	// global MPI communication settings
	extern int grank;
	extern int gsize;
	extern MPI_Comm gcomm;

	void set();
	bool divide(int meshDuplication);
	void finish();
}
}
}



#endif /* SRC_ESINFO_MPIINFO_H_ */
