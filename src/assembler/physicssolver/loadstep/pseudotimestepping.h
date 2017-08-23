
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_PSEUDOTIMESTEPPING_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_PSEUDOTIMESTEPPING_H_

#include "loadstepsolver.h"

namespace espreso {

class NonLinearSolverConfiguration;

class PseudoTimeStepping: public LoadStepSolver {

public:
	PseudoTimeStepping(TimeStepSolver &timeStepSolver, const NonLinearSolverConfiguration &configuration, double duration);

	Matrices updateStructuralMatrices(Step &step, Matrices matrices);
	Matrices reassembleStructuralMatrices(Step &step, Matrices matrices);

protected:
	void runNextTimeStep(Step &step);
	void processTimeStep(Step &step);

	const NonLinearSolverConfiguration &_configuration;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_PSEUDOTIMESTEPPING_H_ */
