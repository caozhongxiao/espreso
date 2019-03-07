
#ifndef SRC_PHYSICS_SOLVER_TIMESTEP_TIMESTEPSOLVER_H_
#define SRC_PHYSICS_SOLVER_TIMESTEP_TIMESTEPSOLVER_H_

#include <string>

namespace espreso {

struct LoadStepSolverConfiguration;
class LoadStepSolver;
class Assembler;

class TimeStepSolver {

	friend class LoadStepSolver;

public:
	TimeStepSolver(Assembler &assembler): _assembler(assembler) {}
	virtual ~TimeStepSolver() {}

	virtual bool hasSameMode(const LoadStepSolverConfiguration &configuration) const =0;
	virtual void solve(LoadStepSolver &loadStepSolver) =0;
	virtual std::string name() =0;

protected:
	Assembler &_assembler;
};

}



#endif /* SRC_PHYSICS_SOLVER_TIMESTEP_TIMESTEPSOLVER_H_ */
