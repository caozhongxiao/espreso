
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREGMRES_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREGMRES_H_

#include "hypreboomeramg.h"
#include "hypreparasalis.h"

namespace espreso {

struct HYPREGMRESConfiguration: public ECFObject {

	enum class PRECONDITIONER {
		BoomerAMG,
		Parasalis,
		NONE
	};
	PRECONDITIONER preconditioner;

	HYPREBoomerAMGConfiguration boomeramg;
	HYPREParasalisConfiguration parasalis;

	double relative_conv_tol, absolute_conv_tol;
	int max_iterations, restarts;
	
	enum class SOLVER_INFO {
		NO_INFO,
		SETUP_INFO,
		SOLVE_INFO,
		SETUP_SOLVE_INFO
	};
	SOLVER_INFO solver_info;
            
	HYPREGMRESConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREGMRES_H_ */
