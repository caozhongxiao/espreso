
#ifndef SRC_WRAPPERS_METIS_WPARMETIS_H_
#define SRC_WRAPPERS_METIS_WPARMETIS_H_

namespace espreso {

struct ParMETIS {

	enum class METHOD {
		ParMETIS_V3_PartKway,
		ParMETIS_V3_RefineKway,
		ParMETIS_V3_AdaptiveRepart
	};

	static eslocal call(
			METHOD method,
			eslocal *edistribution,
			eslocal *eframes, eslocal *eneighbors,
			eslocal dimensions, double *coordinates,
			eslocal verticesWeightCount, eslocal *verticesWeights, eslocal *edgeWeights,
			eslocal *partition);
};

}



#endif /* SRC_WRAPPERS_METIS_WPARMETIS_H_ */
