
#ifndef SRC_LINEARSOLVER_MULTIGRID_MULTIGRID_H_
#define SRC_LINEARSOLVER_MULTIGRID_MULTIGRID_H_

#include "../linearsolver.h"
#include "../../config/ecf/linearsolver/multigrid.h"

namespace espreso {

class Instance;
class HypreData;

class MultigridSolver: public LinearSolver {
public:

	MultigridSolver(Instance *instance, const MultigridConfiguration &configuration);
	virtual ~MultigridSolver();

	void update(Matrices matrices);
	void solve();
	void finalize();

	bool glueDomainsByLagrangeMultipliers() const { return false; }
	bool applyB1Scaling() const {return false;}
	bool applyB1LagrangeRedundancy()const { return false;}

	double& precision();

protected:
	Instance *_instance;
	MultigridConfiguration _configuration;

	HypreData *_hypreData;
};

}

#endif /* SRC_LINEARSOLVER_MULTIGRID_MULTIGRID_H_ */
