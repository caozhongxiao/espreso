
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_

#include "transientsolver.h"
#include "nonlinearsolver.h"
#include "../../solver/feti.h"
#include "../../solver/multigrid.h"

namespace espreso {

struct LoadStepConfiguration: public ECFObject {

	enum class TYPE {
		STEADY_STATE,
		TRANSIENT
	};

	enum class MODE {
		LINEAR,
		NONLINEAR
	};

	enum class SOLVER {
		FETI,
		MULTIGRID
	};

	double duration_time;

	TYPE type;
	MODE mode;
	SOLVER solver;

	NonLinearSolverConfiguration nonlinear_solver;
	TransientSolverConfiguration transient_solver;

	FETISolverConfiguration feti;
	MultigridConfiguration multigrid;

	LoadStepConfiguration(const std::string &firstResidualName, const std::string &secondResidualName);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_ */
