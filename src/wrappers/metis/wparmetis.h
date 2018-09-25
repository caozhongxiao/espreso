
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

	static eslocal call(
			METHOD method,
			MPISubset &subset,
			std::vector<eslocal> &eframes, std::vector<eslocal> &eneighbors,
			std::vector<eslocal> &partition);

	static eslocal call(
			METHOD method,
			MPISubset &subset,
			eslocal *edistribution,
			eslocal *eframes, eslocal *eneighbors,
			eslocal dimensions, double *coordinates,
			eslocal verticesWeightCount, eslocal *verticesWeights, eslocal *edgeWeights,
			eslocal *partition);
};

}



#endif /* SRC_WRAPPERS_METIS_WPARMETIS_H_ */
