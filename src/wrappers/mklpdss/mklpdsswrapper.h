
#ifndef SRC_WRAPPERS_MKLPDSS_MKLPDSSWRAPPER_H_
#define SRC_WRAPPERS_MKLPDSS_MKLPDSSWRAPPER_H_

namespace espreso {

enum class MatrixType : int;
struct MKLPDSSConfiguration;
struct MKLPDSSDataHolder;

class MKLPDSSData {
public:
	MKLPDSSData(esint nrows);

	void insertK(MatrixType mtype, esint *rowPtrs, esint *colIndices, double *values);
	void insertRHS(double *rhsValues);

	void solve(const MKLPDSSConfiguration &configuration, double *solution);

	~MKLPDSSData();

protected:
	void call(esint phase);

	esint _roffset, _nrows;

	MKLPDSSDataHolder *_data;
};

}



#endif /* SRC_WRAPPERS_MKLPDSS_MKLPDSSWRAPPER_H_ */
