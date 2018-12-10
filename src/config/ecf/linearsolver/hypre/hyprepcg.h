
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPCG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPCG_H_

#include "hypreboomeramg.h"

namespace espreso {

struct HYPREPCGConfiguration: public ECFObject {

	enum class PRECONDITIONER {
		BoomerAMG,
		Parasalis
	};

	double relative_conv_tol;

	PRECONDITIONER preconditioner;

	HYPREBoomerAMGConfiguration boomeramg;

	HYPREPCGConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPCG_H_ */
