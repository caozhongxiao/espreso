
#ifndef SRC_ESINFO_ENVINFO_H_
#define SRC_ESINFO_ENVINFO_H_

#include "mpi.h"
#include <cstddef>

namespace espreso {
namespace info {
namespace env {
	extern size_t MKL_NUM_THREADS;
	extern size_t OMP_NUM_THREADS;
	extern size_t SOLVER_NUM_THREADS;
	extern size_t PAR_NUM_THREADS;

	void set();
}
}
}



#endif /* SRC_ESINFO_ENVINFO_H_ */
