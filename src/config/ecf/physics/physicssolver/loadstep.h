
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_

#include "transientfirstorderimplicit.h"
#include "transientsecondorderimplicit.h"
#include "nonlinearsolver.h"
#include "config/ecf/linearsolver/feti.h"
#include "config/ecf/linearsolver/hypre/hypre.h"
#include "config/ecf/linearsolver/mklpdss.h"

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
		HYPRE,
		MKLPDSS
	};

	double duration_time;

	TYPE type;
	MODE mode;
	SOLVER solver;

	NonLinearSolverConfiguration nonlinear_solver;
	TransientFirstOrderImplicitConfiguration transient_solver;
	TransientSecondOrderImplicitConfiguration sm_transient_solver;

	FETISolverConfiguration feti;
	HypreConfiguration hypre;
	MKLPDSSConfiguration mklpdss;

	LoadStepConfiguration(const std::string &firstResidualName, const std::string &secondResidualName);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_ */
