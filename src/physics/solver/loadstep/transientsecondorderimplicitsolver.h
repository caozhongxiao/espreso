
#ifndef SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTSECONDORDERIMPLICITSOLVER_H_
#define SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTSECONDORDERIMPLICITSOLVER_H_

#include "loadstepsolver.h"

namespace espreso {

class NodeData;
class TransientSecondOrderImplicitConfiguration;

class TransientSecondOrderImplicit: public LoadStepSolver {

public:
	TransientSecondOrderImplicit(TransientSecondOrderImplicit *previous, Assembler &assembler, TimeStepSolver &timeStepSolver, TransientSecondOrderImplicitConfiguration &configuration, double duration);

	bool hasSameType(const LoadStepConfiguration &configuration) const;
	std::string name();

	Matrices updateStructuralMatrices(Matrices matrices);

protected:
	void initLoadStep();
	void runNextTimeStep();
	void processTimeStep();

	void updateConstants();

	TransientSecondOrderImplicitConfiguration &_configuration;
	double _alpha, _delta;
	double _massDumping, _stiffnessDumping;
	double _nTimeShift;
	double _newmarkConsts[8];

	NodeData *U, *dU, *V, *W, *X, *Y, *Z, *dTK, *dTM;
};

}



#endif /* SRC_PHYSICS_SOLVER_LOADSTEP_TRANSIENTSECONDORDERIMPLICITSOLVER_H_ */
