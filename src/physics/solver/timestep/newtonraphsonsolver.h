
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
	bool hasSameMode(const LoadStepSolverConfiguration &configuration) const;

	void solve(LoadStepSolver &loadStepSolver);
	std::string name();
	void setSolverParams();

protected:
	static const char *statusNOT, *statusYES;
	NonLinearSolverConfiguration &_configuration;

	NodeData *_solution;

	double _lsAlpha;
	double _lsInc;
	double _lsSInc;
	double _firstConvergenceValue;
	double _firstCriterionValue;
	double _secondConvergenceValue;
	double _secondCriterionValue;
	const char *_status;
};

}



#endif /* SRC_PHYSICS_SOLVER_TIMESTEP_NEWTONRAPHSONSOLVER_H_ */
