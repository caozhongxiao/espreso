
#ifndef SRC_LINEARSOLVER_LINEARSOLVER_H_
#define SRC_LINEARSOLVER_LINEARSOLVER_H_

namespace espreso {

enum Matrices: int;

class LinearSolver {

public:
	void solve(Matrices matrices);

	virtual bool glueDomainsByLagrangeMultipliers() const =0;
	virtual bool applyB1Scaling() const =0;
	virtual bool applyB1LagrangeRedundancy() const =0;

	virtual double& precision() =0;

	virtual ~LinearSolver() {}

protected:
	virtual void update(Matrices matrices) =0;
	virtual void solve() =0;
	virtual void finalize() =0;
};

}



#endif /* SRC_LINEARSOLVER_LINEARSOLVER_H_ */
