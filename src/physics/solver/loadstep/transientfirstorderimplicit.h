
#ifndef SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_
#define SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_

#include "../../solver/loadstep/loadstepsolver.h"

namespace espreso {

class NodeData;
class TransientSolverConfiguration;

class TransientFirstOrderImplicit: public LoadStepSolver {

public:
	TransientFirstOrderImplicit(Assembler &assembler, TimeStepSolver &timeStepSolver, TransientSolverConfiguration &configuration, double duration);

	std::string name();

	Matrices updateStructuralMatrices(Matrices matrices);
	Matrices reassembleStructuralMatrices(Matrices matrices);

protected:
	void initLoadStep();
	void runNextTimeStep();
	void processTimeStep();

	TransientSolverConfiguration &_configuration;
	double _alpha;
	double _nTimeStep;

	static size_t loadStep;

	NodeData *U, *dU, *V, *X, *Y, *dTK;
};

}



#endif /* SRC_ASSEMBLER_PHYSICSSOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICIT_H_ */
