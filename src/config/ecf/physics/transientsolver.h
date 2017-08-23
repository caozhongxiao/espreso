
#ifndef SRC_CONFIG_ECF_PHYSICS_TRANSIENTSOLVER_H_
#define SRC_CONFIG_ECF_PHYSICS_TRANSIENTSOLVER_H_

#include "../../configuration.h"

namespace espreso {

struct AutoTimeSteppingConfiguration: public ECFObject {

	bool allowed;

	double initial_time_step, min_time_step, max_time_step;

	AutoTimeSteppingConfiguration();
};


struct TransientSolverConfiguration: public ECFObject {

	enum class METHOD {
		CRANK_NICOLSON,
		FORWARD_DIFF,
		GALERKIN,
		BACKWARD_DIFF,
		USER
	};

	METHOD method;
	AutoTimeSteppingConfiguration auto_time_stepping;
	double alpha, time_step;

	TransientSolverConfiguration();
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_TRANSIENTSOLVER_H_ */
