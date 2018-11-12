
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_TIMESTEPSOLVER_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_TIMESTEPSOLVER_H_

#include <string>

namespace espreso {

class LinearSolver;
class LoadStepSolver;
class Composer;
struct Step;

class TimeStepSolver {

	friend class LoadStepSolver;

public:
	TimeStepSolver(const std::string &description, Composer &assembler);
	virtual ~TimeStepSolver() {}

	virtual void solve(LoadStepSolver &loadStepSolver) =0;

	std::string description() const;

protected:
	std::string _description;
	Composer &_composer;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_TIMESTEPSOLVER_H_ */
