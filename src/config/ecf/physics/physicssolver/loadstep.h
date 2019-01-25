
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_

#include "transientsolver.h"
#include "nonlinearsolver.h"
#include "config/ecf/linearsolver/feti.h"
#include "config/ecf/linearsolver/hypre/hypre.h"

namespace espreso {

struct LoadStepConfiguration: public ECFDescription {

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
		HYPRE
	};

	double duration_time;

	TYPE type;
	MODE mode;
	SOLVER solver;

	NonLinearSolverConfiguration nonlinear_solver;
	TransientSolverConfiguration transient_solver;

	FETISolverConfiguration feti;
	HypreConfiguration hypre;

	LoadStepConfiguration(const std::string &firstResidualName, const std::string &secondResidualName);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_ */
