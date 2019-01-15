
#ifndef SRC_ESINFO_MPIINFO_H_
#define SRC_ESINFO_MPIINFO_H_

#include "mpi.h"

namespace espreso {
namespace info {
namespace mpi {
	extern int MPIrank;
	extern int MPIsize;
	extern MPI_Comm MPICommunicator;

	void setMPI();
}
}
}



#endif /* SRC_ESINFO_MPIINFO_H_ */
