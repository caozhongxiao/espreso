
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_NEWTONRHAPSON_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_NEWTONRHAPSON_H_

#include "timestepsolver.h"

#include <vector>

namespace espreso {

class NonLinearSolverConfiguration;

class NewtonRhapson: public TimeStepSolver {

public:
	NewtonRhapson(Assembler &assembler, const NonLinearSolverConfiguration &configuration);

	void solve(Step &step, LoadStepSolver &loadStepSolver);

protected:
	const NonLinearSolverConfiguration &_configuration;

	std::vector<std::vector<double> > _solution;
	std::vector<std::vector<double> > _f_ext;
	std::vector<std::vector<double> > _f_R_BtLambda;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_NEWTONRHAPSON_H_ */
