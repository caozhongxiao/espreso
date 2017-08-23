
#ifndef SRC_CONFIG_ECF_PHYSICS_LOADSTEPS_H_
#define SRC_CONFIG_ECF_PHYSICS_LOADSTEPS_H_

#include "transientsolver.h"
#include "nonlinearsolver.h"
#include "../solver/feti.h"
#include "../solver/multigrid.h"

namespace espreso {

struct LoadStepsConfiguration: public ECFObject {

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

	TransientSolverConfiguration transient_solver;

	FETISolverConfiguration feti;
	MultigridConfiguration multigrid;

	LoadStepsConfiguration();
};

struct AdvectionDiffusionLoadStepsConfiguration: public LoadStepsConfiguration {
	AdvectionDiffusionNonLinearSolverConfiguration nonlinear_solver;

	AdvectionDiffusionLoadStepsConfiguration();
};

struct StructuralMechanicsLoadStepsConfiguration: public LoadStepsConfiguration {
	StructuralMechanicsNonLinearSolverConfiguration nonlinear_solver;

	StructuralMechanicsLoadStepsConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_PHYSICS_LOADSTEPS_H_ */
