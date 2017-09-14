
#ifndef SRC_CONFIG_ECF_SOLVER_MULTIGRID_H_
#define SRC_CONFIG_ECF_SOLVER_MULTIGRID_H_

#include "../../configuration.h"

namespace espreso {

enum class MULTIGRID_SOLVER {
	CG = 0,
	GMRES = 1,
	FGMRES = 2,
	BICGS = 3,
	BICGSTAB = 4,
	TFQMR = 5,
	SYMQMR = 6,
	SUPERLU = 7,
	SUPERLUX = 8
};

enum class MULTIGRID_PRECONDITIONER {
	DIAGONAL = 0,
	PILUT = 1,
	EUCLID = 2,
	PARASAILS = 3,
	BOOMERAMG = 4,
	POLY = 5,
	MLI = 6
};

struct MultigridConfiguration: public ECFObject {

	double precision;
	size_t max_iterations;

	MULTIGRID_SOLVER solver;
	MULTIGRID_PRECONDITIONER preconditioner;

	MultigridConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_SOLVER_MULTIGRID_H_ */
