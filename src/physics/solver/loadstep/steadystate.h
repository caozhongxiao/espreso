
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_STEADYSTATE_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_STEADYSTATE_H_

#include "loadstepsolver.h"

namespace espreso {

class SteadyStateSolver: public LoadStepSolver {

public:
	SteadyStateSolver(Assembler &assembler, TimeStepSolver &timeStepSolver, double duration);

	std::string name();

	Matrices updateStructuralMatrices(Matrices matrices);

protected:
	void runNextTimeStep();
	void processTimeStep();
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_STEADYSTATE_H_ */
