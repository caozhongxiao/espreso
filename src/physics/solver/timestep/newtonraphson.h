
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_NEWTONRAPHSON_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_NEWTONRAPHSON_H_

#include <vector>
#include "../../solver/timestep/timestepsolver.h"

namespace espreso {

class NonLinearSolverConfiguration;

class NewtonRaphson: public TimeStepSolver {

public:
	NewtonRaphson(Assembler &assembler, LinearSolver &solver, NonLinearSolverConfiguration &configuration);

	void solve(LoadStepSolver &loadStepSolver);
	std::string name();

protected:
	NonLinearSolverConfiguration &_configuration;

	std::vector<std::vector<double> > _solution;
	std::vector<std::vector<double> > _f_ext;
	std::vector<std::vector<double> > _f_R_BtLambda;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_NEWTONRAPHSON_H_ */
