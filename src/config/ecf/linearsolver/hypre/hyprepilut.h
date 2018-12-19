
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPILUT_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPILUT_H_

#include "../../../configuration.h"

namespace espreso {

struct HYPREPilutConfiguration: public ECFObject {

    int max_iter;
    double drop_tol;
    int row_size;
	
	HYPREPilutConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREPILUT_H_ */
