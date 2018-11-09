
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_STEADYSTATE_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_STEADYSTATE_H_

#include "../../solver/loadstep/loadstepsolver.h"

namespace espreso {

class SteadyStateSolver: public LoadStepSolver {

public:
	SteadyStateSolver(TimeStepSolver &timeStepSolver, double duration);

	Matrices updateStructuralMatrices(Matrices matrices);
	Matrices reassembleStructuralMatrices(Matrices matrices);

protected:
	void runNextTimeStep();
	void processTimeStep();
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_STEADYSTATE_H_ */
