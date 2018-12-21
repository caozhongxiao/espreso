
#ifndef SRC_WRAPPERS_METIS_WPARMETIS_H_
#define SRC_WRAPPERS_METIS_WPARMETIS_H_

#include <vector>

namespace espreso {

struct MPISubset;

struct ParMETIS {

	enum class METHOD {
		ParMETIS_V3_PartKway,
		ParMETIS_V3_RefineKway,
		ParMETIS_V3_AdaptiveRepart,
		ParMETIS_V3_PartGeom
	};

	static esint call(
			METHOD method,
			MPISubset &subset,
			std::vector<esint> &eframes, std::vector<esint> &eneighbors,
			std::vector<esint> &partition);

	static esint call(
			METHOD method,
			MPISubset &subset,
			esint *edistribution,
			esint *eframes, esint *eneighbors,
			esint dimensions, double *coordinates,
			esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
			esint *partition);
};

}



#endif /* SRC_WRAPPERS_METIS_WPARMETIS_H_ */
