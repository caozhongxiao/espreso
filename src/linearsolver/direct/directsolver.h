
#ifndef SRC_LINEARSOLVER_DIRECT_DIRECTSOLVER_H_
#define SRC_LINEARSOLVER_DIRECT_DIRECTSOLVER_H_

#include "linearsolver/linearsolver.h"

namespace espreso {

class DirectSolver: public LinearSolver {
public:

	DirectSolver(DataHolder *data);
	virtual ~DirectSolver();

	void update(Matrices matrices);
	void solve();
	void finalize();

	double& precision();
};

}



#endif /* SRC_LINEARSOLVER_DIRECT_DIRECTSOLVER_H_ */
