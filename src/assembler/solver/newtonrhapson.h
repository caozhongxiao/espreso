
#ifndef SRC_ASSEMBLER_SOLVER_NEWTONRHAPSON_H_
#define SRC_ASSEMBLER_SOLVER_NEWTONRHAPSON_H_

#include "solver.h"

namespace espreso {

struct NonLinearSolverBase;

class NewtonRhapson: public Solver
{
public:
	NewtonRhapson(
			Mesh *mesh,
			Physics* physics,
			LinearSolver* linearSolver,
			store::ResultStore* store,
			const NonLinearSolverBase &configuration,
			Matrices restriction = Matrices::NONE);

	virtual void run(Step &step);

	virtual void init(Step &step);
	virtual void preprocess(Step &step);
	virtual void solve(Step &step);
	virtual void postprocess(Step &step);
	virtual void finalize(Step &step);

protected:
	const NonLinearSolverBase &_configuration;
};

}



#endif /* SRC_ASSEMBLER_SOLVER_NEWTONRHAPSON_H_ */
