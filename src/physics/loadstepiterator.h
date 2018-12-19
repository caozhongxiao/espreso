
#ifndef SRC_PHYSICS_LOADSTEPITERATOR_H_
#define SRC_PHYSICS_LOADSTEPITERATOR_H_

namespace espreso {

class LoadStepSolver;
class TimeStepSolver;
class Assembler;
class LinearSolver;

class LoadStepIterator {

public:
	LoadStepIterator();

	bool next();

protected:

	template <typename TPhysics> bool next(TPhysics &configuration);

	LoadStepSolver *_loadStepSolver;
	TimeStepSolver *_timeStepSolver;
	Assembler *_assembler;
	LinearSolver * _linearSolver;
};

}



#endif /* SRC_PHYSICS_LOADSTEPITERATOR_H_ */
