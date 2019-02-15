
#ifndef SRC_WRAPPERS_HYPRE_HYPREWRAPPER_H_
#define SRC_WRAPPERS_HYPRE_HYPREWRAPPER_H_

namespace espreso {

struct HypreConfiguration;
struct HYPREDataHolder;

class HypreData {
public:
	HypreData(esint nrows);

	void insertCSR(esint nrows, esint offset, esint *rowPrts, esint *colIndices, double *values, double *rhsValues);
	void finalizePattern();

	void solve(const HypreConfiguration &configuration, esint nrows, double *solution);

	~HypreData();

protected:
	esint _roffset, _nrows;

	HYPREDataHolder *_data;

	bool _finalized;
};

}

#endif /* SRC_WRAPPERS_HYPRE_HYPREWRAPPER_H_ */
