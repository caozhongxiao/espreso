
#ifndef SRC_PHYSICS_SOLVER_TIMESTEP_LINEAR_H_
#define SRC_PHYSICS_SOLVER_TIMESTEP_LINEAR_H_

#include "timestepsolver.h"

namespace espreso {

class LinearTimeStep: public TimeStepSolver {

public:
	LinearTimeStep(Assembler &assembler);

	void solve(LoadStepSolver &loadStepSolver);
	std::string name();

};

}

#endif /* SRC_PHYSICS_SOLVER_TIMESTEP_LINEAR_H_ */
