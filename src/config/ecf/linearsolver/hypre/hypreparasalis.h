
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPARASALIS_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPARASALIS_H_

#include "../../../configuration.h"

namespace espreso {

struct HYPREParasalisConfiguration: public ECFObject {

	double convergence_tolerance;

	HYPREParasalisConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPARASALIS_H_ */
