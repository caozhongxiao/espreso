
#ifndef SRC_LINEARSOLVER_MULTIGRID_MULTIGRID_H_
#define SRC_LINEARSOLVER_MULTIGRID_MULTIGRID_H_

#include "../linearsolver.h"
#include "../../assembler/instance.h"
#include "../../basis/logging/timeeval.h"
#include "../../config/ecf/solver/multigrid.h"

namespace espreso {

class MultigridSolver : public LinearSolver {
public:

	MultigridSolver(Instance *instance, const MultigridConfiguration &configuration);
	virtual ~MultigridSolver();

	void update(Matrices matrices);
	void solve();
	void finalize();

	bool glueDomainsByLagrangeMultipliers() const { return false; }
	bool applyB1Scaling() const {return false;}
	bool applyB1LagrangeRedundancy()const { return false;}

	double& precision() { return configuration.precision; }

private:
	Instance *instance;
	MultigridConfiguration configuration;
	TimeEval timeEvalMain;
};

}

#endif /* SRC_LINEARSOLVER_MULTIGRID_MULTIGRID_H_ */
