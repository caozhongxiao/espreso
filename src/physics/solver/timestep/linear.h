
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_LINEAR_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_LINEAR_H_

#include "../../solver/timestep/timestepsolver.h"

namespace espreso {

class LinearTimeStep: public TimeStepSolver {

public:
	LinearTimeStep(Composer &composer);

	void solve(LoadStepSolver &loadStepSolver);

};

}

#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_LINEAR_H_ */
