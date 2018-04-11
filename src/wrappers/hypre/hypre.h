/*
 * hypre.h
 *
 *  Created on: Apr 9, 2018
 *      Author: beh01
 */

#ifndef SRC_WRAPPERS_HYPRE_HYPRE_H_
#define SRC_WRAPPERS_HYPRE_HYPRE_H_

#include "mpi.h"

namespace espreso {

struct  HYPRE {
public:

	static void Solve(int rank, int processes, MPI_Comm communicator);

};

}
#endif /* SRC_WRAPPERS_HYPRE_HYPRE_H_ */
