
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_LOADSTEPSOLVER_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_LOADSTEPSOLVER_H_

#include <string>

namespace espreso {

struct HeatTransferConfiguration;
struct StructuralMechanicsConfiguration;
class TimeStepSolver;
class Provider;
enum Matrices: int;

class LoadStepSolver {

	friend class LinearTimeStep;
	friend class NewtonRaphson;

public:
//	static LoadStepSolver* create(const HeatTransferConfiguration &configuration);
//	static LoadStepSolver* create(const StructuralMechanicsConfiguration &configuration);

	LoadStepSolver(const std::string &description, TimeStepSolver &timeStepSolver, double duration);
	virtual ~LoadStepSolver() {}

	void run();

	std::string description() const;
	double duration() const;

protected:
	virtual Matrices updateStructuralMatrices(Matrices matrices) =0;
	virtual Matrices reassembleStructuralMatrices(Matrices matrices) =0;

	virtual void initLoadStep();
	virtual bool hasNextTimeStep();
	virtual void runNextTimeStep() =0;
	virtual void processTimeStep() =0;
	virtual void finalizeLoadStep();

	TimeStepSolver *_timeStep;

	std::string _description;
	TimeStepSolver &_timeStepSolver;
	Provider &_composer;
	double _duration;

	double _startTime;
	double _precision;
};

}


#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_LOADSTEPSOLVER_H_ */
