
#ifndef SRC_LINEARSOLVER_MKLPDSS_MKLPDSSSOLVER_H_
#define SRC_LINEARSOLVER_MKLPDSS_MKLPDSSSOLVER_H_

#include "linearsolver/linearsolver.h"

namespace espreso {

struct MKLPDSSConfiguration;
struct MKLPDSSData;

class MKLPDSSSolver: public LinearSolver {
public:

	MKLPDSSSolver(DataHolder *data, MKLPDSSConfiguration &configuration);
	virtual ~MKLPDSSSolver();

	void update(Matrices matrices);
	void solve();
	void finalize();

	double& precision() { return _precision; }

protected:
	MKLPDSSConfiguration &_configuration;

	MKLPDSSData *_mklpdssData;

	double _precision; // dummy since MKL PDSS does not have this option
};

}



#endif /* SRC_LINEARSOLVER_MKLPDSS_MKLPDSSSOLVER_H_ */
