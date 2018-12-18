
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREFlexGMRES_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREFlexGMRES_H_

#include "hypreboomeramg.h"
#include "hypreparasails.h"
#include "hypreeuclid.h"

namespace espreso {

struct HYPREFlexGMRESConfiguration: public ECFObject {

	enum class PRECONDITIONER {
		BoomerAMG,
		ParaSails,
		Euclid,
		NONE
	};
	PRECONDITIONER preconditioner;

	HYPREBoomerAMGConfiguration boomeramg;
	HYPREParaSailsConfiguration parasails;
	HYPREEuclidConfiguration euclid;

	double relative_conv_tol, absolute_conv_tol;
	int max_iterations, restarts;
	
	enum class SOLVER_INFO {
		NO_INFO,
		SETUP_INFO,
		SOLVE_INFO,
		SETUP_SOLVE_INFO
	};
	SOLVER_INFO solver_info;
            
	HYPREFlexGMRESConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREFlexGMRES_H_ */
