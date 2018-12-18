
#ifndef SRC_PHYSICS_SOLVER_TIMESTEP_TIMESTEPSOLVER_H_
#define SRC_PHYSICS_SOLVER_TIMESTEP_TIMESTEPSOLVER_H_

#include <string>

namespace espreso {

class LoadStepSolver;
class Assembler;
class LinearSolver;

class TimeStepSolver {

	friend class LoadStepSolver;

public:
	TimeStepSolver(Assembler &assembler, LinearSolver &solver): _assembler(assembler), _solver(solver) {}
	virtual ~TimeStepSolver() {}

	virtual void solve(LoadStepSolver &loadStepSolver) =0;
	virtual std::string name() =0;

protected:
	Assembler &_assembler;
	LinearSolver &_solver;
};

}



#endif /* SRC_PHYSICS_SOLVER_TIMESTEP_TIMESTEPSOLVER_H_ */
