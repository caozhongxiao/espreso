
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPARASALIS_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPARASALIS_H_

#include "config/description.h"

namespace espreso {

struct HYPREParasalisConfiguration: public ECFDescription {

	double convergence_tolerance;

	HYPREParasalisConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPARASALIS_H_ */
