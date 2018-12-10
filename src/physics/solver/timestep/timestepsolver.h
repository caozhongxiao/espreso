
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_TIMESTEPSOLVER_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_TIMESTEPSOLVER_H_

#include <string>

namespace espreso {

class LinearSolver;
class LoadStepSolver;
class Provider;

class TimeStepSolver {

	friend class LoadStepSolver;

public:
//	static TimeStepSolver* create(const HeatTransferConfiguration &configuration);
//	static TimeStepSolver* create(const StructuralMechanicsConfiguration &configuration);

	TimeStepSolver(const std::string &description, Provider &assembler);
	virtual ~TimeStepSolver() {}

	virtual void solve(LoadStepSolver &loadStepSolver) =0;

	std::string description() const;

protected:
	std::string _description;
	Provider &_composer;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_TIMESTEP_TIMESTEPSOLVER_H_ */
