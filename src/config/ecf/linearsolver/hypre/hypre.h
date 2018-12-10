
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPRE_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPRE_H_

#include "hypreboomeramg.h"
#include "hyprepcg.h"

namespace espreso {

struct HypreConfiguration: public ECFObject {

	enum class SOLVER_TYPE {
		BoomerAMG,
		PCG,
		GMRES,
	};

	SOLVER_TYPE solver_type;

	HYPREBoomerAMGConfiguration boomeramg;
	HYPREPCGConfiguration pcg;

	HypreConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPRE_H_ */
