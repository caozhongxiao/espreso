
#ifndef SRC_PHYSICS_SOLVER_TIMESTEP_LINEARTIMESOLVER_H_
#define SRC_PHYSICS_SOLVER_TIMESTEP_LINEARTIMESOLVER_H_

#include "timestepsolver.h"

namespace espreso {

class LinearTimeStep: public TimeStepSolver {

public:
	LinearTimeStep(LinearTimeStep *previous, Assembler &assembler);
	bool hasSameMode(const LoadStepConfiguration &configuration) const;

	void solve(LoadStepSolver &loadStepSolver);
	std::string name();

};

}

#endif /* SRC_PHYSICS_SOLVER_TIMESTEP_LINEARTIMESOLVER_H_ */
