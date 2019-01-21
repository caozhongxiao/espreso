
#ifndef SRC_ESINFO_MPIINFO_H_
#define SRC_ESINFO_MPIINFO_H_

#include "mpi.h"

namespace espreso {
namespace info {
namespace mpi {
	extern int rank;
	extern int size;
	extern MPI_Comm comm;

	void set();
}
}
}



#endif /* SRC_ESINFO_MPIINFO_H_ */
