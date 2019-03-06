
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_STEADYSTATESOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_STEADYSTATESOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

class SteadyStateSolver: public LoadStepSolver {

public:
	SteadyStateSolver(SteadyStateSolver *previous, Assembler &assembler, TimeStepSolver &timeStepSolver, double duration);

	bool hasSameType(const LoadStepConfiguration &configuration) const;
	std::string name();

	Matrices updateStructuralMatrices(Matrices matrices);

protected:
	void runNextTimeStep();
	void processTimeStep();
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_STEADYSTATESOLVER_H_ */
