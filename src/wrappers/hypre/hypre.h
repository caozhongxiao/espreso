
#ifndef SRC_WRAPPERS_HYPRE_HYPRE_H_
#define SRC_WRAPPERS_HYPRE_HYPRE_H_

#include "mpi.h"

#include "include/HYPRE_IJ_mv.h"

namespace espreso {

class SparseMatrix;
struct MultigridConfiguration;

class HypreData {
	friend class HYPRE;
public:
	HypreData(MPI_Comm &comm, eslocal nrows);

	void insertCSR(eslocal nrows, eslocal offset, eslocal *rowPrts, eslocal *colIndices, double *values, double *rhsValues);
	void insertIJV(eslocal nrows, eslocal offset, eslocal size, eslocal *rowIndices, eslocal *colIndices, double *values, double *rhsValues);
	void finalizePattern();

	~HypreData();

protected:
	MPI_Comm &_comm;
	eslocal _roffset, _nrows;

	HYPRE_IJMatrix _K;
	HYPRE_IJVector _f, _x;

	bool _finalized;
};

struct HYPRE {
	static void solve(const MultigridConfiguration &configuration, HypreData &data, eslocal nrows, double *solution);
};

}
#endif /* SRC_WRAPPERS_HYPRE_HYPRE_H_ */
