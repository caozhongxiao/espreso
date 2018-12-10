
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREBOOMERAMG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREBOOMERAMG_H_

#include "../../../configuration.h"

namespace espreso {

struct HYPREBoomerAMGConfiguration: public ECFObject {

	enum class CYCLE_TYPE {
		V_CYCLE,
		W_CYCLE
	};

	double convergence_tolerance;
	int min_iterations, max_iterations;

	CYCLE_TYPE cycle_type;

	HYPREBoomerAMGConfiguration();
};


}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREBOOMERAMG_H_ */
