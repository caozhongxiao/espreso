
#ifndef SRC_LINEARSOLVER_LINEARSOLVER_H_
#define SRC_LINEARSOLVER_LINEARSOLVER_H_

namespace espreso {

struct DataHolder;
enum Matrices: int;

class LinearSolver {

public:
	LinearSolver(DataHolder *data);

	void solve(Matrices matrices);

	virtual double& precision() =0;

	virtual ~LinearSolver() {}

protected:
	virtual void update(Matrices matrices) =0;
	virtual void solve() =0;
	virtual void finalize() =0;

	DataHolder *_data;
};

}



#endif /* SRC_LINEARSOLVER_LINEARSOLVER_H_ */
