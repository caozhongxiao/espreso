
#ifndef SRC_WRAPPERS_HYPRE_HYPREWRAPPER_H_
#define SRC_WRAPPERS_HYPRE_HYPREWRAPPER_H_

namespace espreso {

struct HypreConfiguration;
struct HYPREData;

class HypreData {
	friend class HYPRE;
public:
	HypreData(esint nrows);

	void insertCSR(esint nrows, esint offset, esint *rowPrts, esint *colIndices, double *values, double *rhsValues);
	void insertIJV(esint nrows, esint offset, esint size, esint *rowIndices, esint *colIndices, double *values, double *rhsValues);
	void finalizePattern();

	~HypreData();

protected:
	esint _roffset, _nrows;

	HYPREData *_data;

	bool _finalized;
};

struct HYPRE {
	static void solve(const HypreConfiguration &configuration, HypreData &data, esint nrows, double *solution);
};

}
#endif /* SRC_WRAPPERS_HYPRE_HYPREWRAPPER_H_ */
