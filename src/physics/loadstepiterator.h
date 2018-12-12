
#ifndef SRC_PHYSICS_LOADSTEPITERATOR_H_
#define SRC_PHYSICS_LOADSTEPITERATOR_H_

namespace espreso {

struct HeatTransferConfiguration;
struct StructuralMechanicsConfiguration;

class LoadStepSolver;
class TimeStepSolver;
class Assembler;
class LinearSolver;

class LoadStepIterator {

public:
	LoadStepIterator();

	bool next();

protected:

	bool next(HeatTransferConfiguration &configuration);
	bool next(StructuralMechanicsConfiguration &configuration);

	LoadStepSolver *_loadStepSolver;
	TimeStepSolver *_timeStepSolver;
	Assembler *_assembler;
	LinearSolver * _linearSolver;
};

}



#endif /* SRC_PHYSICS_LOADSTEPITERATOR_H_ */
