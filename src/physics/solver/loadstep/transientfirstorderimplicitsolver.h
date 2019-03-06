
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICITSOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICITSOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

class NodeData;
class TransientFirstOrderImplicitConfiguration;

class TransientFirstOrderImplicit: public LoadStepSolver {

public:
	TransientFirstOrderImplicit(TransientFirstOrderImplicit *previous, Assembler &assembler, TimeStepSolver &timeStepSolver, TransientFirstOrderImplicitConfiguration &configuration, double duration);

	bool hasSameType(const LoadStepConfiguration &configuration) const;
	std::string name();

	Matrices updateStructuralMatrices(Matrices matrices);

protected:
	void initLoadStep();
	void runNextTimeStep();
	void processTimeStep();

	TransientFirstOrderImplicitConfiguration &_configuration;
	double _alpha;
	double _nTimeShift;

	NodeData *U, *dU, *V, *X, *Y, *dTK, *dTM;
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTFIRSTORDERIMPLICITSOLVER_H_ */
