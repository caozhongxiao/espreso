
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICITSOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICITSOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

class NodeData;
class TransientFirstOrderImplicitSolverConfiguration;

class TransientFirstOrderImplicit: public LoadStepSolver {

public:
	TransientFirstOrderImplicit(TransientFirstOrderImplicit *previous, Assembler &assembler, TimeStepSolver &timeStepSolver, TransientFirstOrderImplicitSolverConfiguration &configuration, double duration);

	bool hasSameType(const LoadStepSolverConfiguration &configuration) const;
	std::string name();
	void setSolverParams();

	Matrices updateStructuralMatrices(Matrices matrices);

protected:
	static const char *statusINCR, *statusDECR, *statusUNCH;
	void initLoadStep();
	void runNextTimeStep();
	void processTimeStep();

	TransientFirstOrderImplicitSolverConfiguration &_configuration;
	double _alpha;
	double _nTimeShift;

	NodeData *U, *dU, *V, *X, *Y, *dTK, *dTM;

	double _resFreq, _oscilationLimit;
	const char *_status;

};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICITSOLVER_H_ */
