
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_

#include "nonlinear.h"
#include "transientfirstorderimplicit.h"
#include "transientsecondorderimplicit.h"
#include "config/ecf/linearsolver/feti.h"
#include "config/ecf/linearsolver/hypre/hypre.h"
#include "config/ecf/linearsolver/mklpdss.h"

namespace espreso {

struct LoadStepSolverConfiguration: public ECFDescription {

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

	FETISolverConfiguration feti;
	HypreConfiguration hypre;
	MKLPDSSConfiguration mklpdss;

	LoadStepSolverConfiguration();
};

struct HeatTransferLoadStepSolverConfiguration: public LoadStepSolverConfiguration {

	NonLinearSolverConfiguration nonlinear_solver;
	TransientFirstOrderImplicitSolverConfiguration transient_solver;

	HeatTransferLoadStepSolverConfiguration();
};

struct StructuralMechanicsLoadStepSolverConfiguration: public LoadStepSolverConfiguration {

	NonLinearSolverConfiguration nonlinear_solver;
	TransientSecondOrderImplicitSolverConfiguration transient_solver;

	StructuralMechanicsLoadStepSolverConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_LOADSTEP_H_ */
