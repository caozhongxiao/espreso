
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_PSEUDOTIMESTEPPING_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_PSEUDOTIMESTEPPING_H_

#include "loadstepsolver.h"

namespace espreso {

class NonLinearSolverConfiguration;

class PseudoTimeStepping: public LoadStepSolver {

public:
	PseudoTimeStepping(PseudoTimeStepping *previous, Assembler &assembler, TimeStepSolver &timeStepSolver, NonLinearSolverConfiguration &configuration, double duration);

	bool hasSameType(const LoadStepConfiguration &configuration) const;
	std::string name();

	Matrices updateStructuralMatrices(Matrices matrices);

protected:
	void runNextTimeStep();
	void processTimeStep();

	NonLinearSolverConfiguration &_configuration;
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_PSEUDOTIMESTEPPING_H_ */
