
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_

#include "loadstepsolver.h"

namespace espreso {

class NodeData;
class TransientSolverConfiguration;

class TransientFirstOrderImplicit: public LoadStepSolver {

public:
	TransientFirstOrderImplicit(Assembler &assembler, TimeStepSolver &timeStepSolver, TransientSolverConfiguration &configuration, double duration);

	std::string name();

	Matrices updateStructuralMatrices(Matrices matrices);

protected:
	void initLoadStep();
	void runNextTimeStep();
	void processTimeStep();

	TransientSolverConfiguration &_configuration;
	double _alpha;
	double _nTimeStep;

	NodeData *U, *dU, *V, *X, *Y, *dTK;
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_ */
