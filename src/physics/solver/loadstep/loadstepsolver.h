
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_LOADSTEPSOLVER_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_LOADSTEPSOLVER_H_

#include <string>

namespace espreso {

class TimeStepSolver;
class Assembler;
enum Matrices: int;

class LoadStepSolver {

public:

	LoadStepSolver(Assembler &assembler, TimeStepSolver &timeStepSolver, double duration);
	virtual ~LoadStepSolver() {}

	void run();

	virtual std::string name() =0;
	double duration() const;

	virtual Matrices updateStructuralMatrices(Matrices matrices) =0;
	virtual Matrices reassembleStructuralMatrices(Matrices matrices) =0;

protected:
	virtual void initLoadStep();
	virtual bool hasNextTimeStep();
	virtual void runNextTimeStep() =0;
	virtual void processTimeStep() =0;
	virtual void finalizeLoadStep();

	Assembler &_assembler;
	TimeStepSolver &_timeStepSolver;

	double _duration;

	double _startTime;
	double _precision;
};

}


#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_LOADSTEPSOLVER_H_ */