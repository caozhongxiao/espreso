
#ifndef SRC_LINEARSOLVER_HYPRE_HYPRESOLVER_H_
#define SRC_LINEARSOLVER_HYPRE_HYPRESOLVER_H_

#include "linearsolver/linearsolver.h"

namespace espreso {

struct HypreConfiguration;
struct HypreData;

class HYPRESolver: public LinearSolver {
public:

	HYPRESolver(DataHolder *data, HypreConfiguration &configuration);
	virtual ~HYPRESolver();

	void update(Matrices matrices);
	void solve();
	void finalize();

	bool glueDomainsByLagrangeMultipliers() const { return false; }
	bool applyB1Scaling() const {return false;}
	bool applyB1LagrangeRedundancy()const { return false;}

	double& precision();

protected:
	HypreConfiguration &_configuration;

	HypreData *_hypreData;
};

}

#endif /* SRC_LINEARSOLVER_HYPRE_HYPRESOLVER_H_ */
