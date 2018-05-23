/*
 * hypre.h
 *
 *  Created on: Apr 9, 2018
 *      Author: beh01
 */

#ifndef SRC_WRAPPERS_HYPRE_HYPRE_H_
#define SRC_WRAPPERS_HYPRE_HYPRE_H_

#include "mpi.h"

#include "../../config/ecf/solver/multigrid.h"

namespace espreso {

class SparseMatrix;

class HypreRegion {
public :
	/* number of rows to be changed */
	int nrows;
	/* number of elements to be set for each row */
	std::vector<int> ncols;
	/* indices of rows in which values are to be set */
	std::vector<int> rows;
	/* column indices of the values to be set */
	const int* cols;
	/* values to be set */
	const double* values;

	HypreRegion() {

	}
};

struct  HYPRE {

public:

	static void Solve(MultigridConfiguration &configuration,
			int rank, int processes, MPI_Comm communicator,
			esglobal start_row, eslocal num_rows,
			std::vector<HypreRegion> values,
			double* f_data,
			double* result);
};

}
#endif /* SRC_WRAPPERS_HYPRE_HYPRE_H_ */
