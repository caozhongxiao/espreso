
#ifndef SRC_PHYSICS_SOLVER_TIMESTEP_NEWTONRAPHSONSOLVER_H_
#define SRC_PHYSICS_SOLVER_TIMESTEP_NEWTONRAPHSONSOLVER_H_

#include "timestepsolver.h"
#include <vector>

namespace espreso {

class NonLinearSolverConfiguration;
struct NodeData;

class NewtonRaphson: public TimeStepSolver {

public:
	NewtonRaphson(NewtonRaphson *previous, Assembler &assembler, NonLinearSolverConfiguration &configuration);
	bool hasSameMode(const LoadStepConfiguration &configuration) const;

	void solve(LoadStepSolver &loadStepSolver);
	std::string name();

protected:
	NonLinearSolverConfiguration &_configuration;

	NodeData *_solution;
};

}



#endif /* SRC_PHYSICS_SOLVER_TIMESTEP_NEWTONRAPHSONSOLVER_H_ */
