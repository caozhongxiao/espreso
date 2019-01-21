
#ifndef SRC_PHYSICS_SOLVER_TIMESTEP_NEWTONRAPHSON_H_
#define SRC_PHYSICS_SOLVER_TIMESTEP_NEWTONRAPHSON_H_

#include "timestepsolver.h"
#include <vector>

namespace espreso {

class NonLinearSolverConfiguration;
struct NodeData;

class NewtonRaphson: public TimeStepSolver {

public:
	NewtonRaphson(Assembler &assembler, NonLinearSolverConfiguration &configuration);

	void solve(LoadStepSolver &loadStepSolver);
	std::string name();

protected:
	NonLinearSolverConfiguration &_configuration;

	NodeData *_solution, *_RHS;
};

}



#endif /* SRC_PHYSICS_SOLVER_TIMESTEP_NEWTONRAPHSON_H_ */
