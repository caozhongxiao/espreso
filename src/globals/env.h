
#ifndef SRC_GLOBALS_ENV_H_
#define SRC_GLOBALS_ENV_H_

#include "mpi.h"
#include <cstddef>

namespace espreso {

struct env {
	static int MPIrank;
	static int MPIsize;
	static MPI_Comm MPICommunicator;

	static size_t MKL_NUM_THREADS;
	static size_t OMP_NUM_THREADS;
	static size_t SOLVER_NUM_THREADS;
	static size_t PAR_NUM_THREADS;

	static size_t verbose_level;
	static size_t measure_level;

	static void setMPI();
};

}



#endif /* SRC_GLOBALS_ENV_H_ */
