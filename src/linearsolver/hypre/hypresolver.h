
#ifndef SRC_LINEARSOLVER_HYPRE_HYPRESOLVER_H_
#define SRC_LINEARSOLVER_HYPRE_HYPRESOLVER_H_

#include "../linearsolver.h"
#include "../../config/ecf/linearsolver/hypre/hypre.h"

namespace espreso {

class Instance;
class HypreData;

class HYPRESolver: public LinearSolver {
public:

	HYPRESolver(Instance *instance, const HypreConfiguration &configuration);
	virtual ~HYPRESolver();

	void update(Matrices matrices);
	void solve();
	void finalize();

	bool glueDomainsByLagrangeMultipliers() const { return false; }
	bool applyB1Scaling() const {return false;}
	bool applyB1LagrangeRedundancy()const { return false;}

	double& precision();

protected:
	Instance *_instance;
	HypreConfiguration _configuration;

	HypreData *_hypreData;
};

}

#endif /* SRC_LINEARSOLVER_HYPRE_HYPRESOLVER_H_ */
