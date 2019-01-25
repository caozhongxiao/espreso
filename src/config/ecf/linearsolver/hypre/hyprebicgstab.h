
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREBiCGSTAB_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREBiCGSTAB_H_

#include "hypreboomeramg.h"
#include "hypreparasalis.h"

namespace espreso {

struct HYPREBiCGSTABConfiguration: public ECFDescription {

	enum class PRECONDITIONER {
		BoomerAMG,
		Parasalis,
		NONE
	};
	PRECONDITIONER preconditioner;

	HYPREBoomerAMGConfiguration boomeramg;
	HYPREParasalisConfiguration parasalis;

	double relative_conv_tol, absolute_conv_tol;
	int max_iterations;
	
	enum class SOLVER_INFO {
		NO_INFO,
		SETUP_INFO,
		SOLVE_INFO,
		SETUP_SOLVE_INFO
	};
	SOLVER_INFO solver_info;

	HYPREBiCGSTABConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREBiCGSTAB_H_ */
