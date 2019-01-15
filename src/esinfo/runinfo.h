
#ifndef SRC_ESINFO_RUNINFO_H_
#define SRC_ESINFO_RUNINFO_H_

#include "mpi.h"
#include <cstddef>

namespace espreso {
namespace info {
namespace run {
	extern size_t MKL_NUM_THREADS;
	extern size_t OMP_NUM_THREADS;
	extern size_t SOLVER_NUM_THREADS;
	extern size_t PAR_NUM_THREADS;

	extern size_t verbose_level;
	extern size_t measure_level;
}
}
}



#endif /* SRC_ESINFO_RUNINFO_H_ */
