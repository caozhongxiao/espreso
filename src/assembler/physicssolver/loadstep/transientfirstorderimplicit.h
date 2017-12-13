
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_

#include "loadstepsolver.h"

namespace espreso {

class NodeData;
class TransientSolverConfiguration;

class TransientFirstOrderImplicit: public LoadStepSolver {

public:
	TransientFirstOrderImplicit(TimeStepSolver &timeStepSolver, const TransientSolverConfiguration &configuration, double duration);

	Matrices updateStructuralMatrices(Step &step, Matrices matrices);
	Matrices reassembleStructuralMatrices(Step &step, Matrices matrices);

protected:
	void initLoadStep(Step &step);
	void runNextTimeStep(Step &step);
	void processTimeStep(Step &step);

	const TransientSolverConfiguration &_configuration;
	double _alpha;
	double _nTimeStep;

	static size_t loadStep;

	NodeData *U, *dU, *V, *X, *Y, *dTK;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_ */
