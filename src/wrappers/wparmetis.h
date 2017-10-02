
#ifndef SRC_WRAPPERS_WPARMETIS_H_
#define SRC_WRAPPERS_WPARMETIS_H_

namespace espreso {

struct ParMETIS {

	enum class METHOD {
		ParMETIS_V3_PartKway,
		ParMETIS_V3_RefineKway,
		ParMETIS_V3_AdaptiveRepart
	};

	static esglobal call(
			METHOD method,
			esglobal *edistribution,
			esglobal *eframes, esglobal *eneighbors,
			esglobal dimensions, double *coordinates,
			esglobal verticesWeightCount, esglobal *verticesWeights, esglobal *edgeWeights,
			esglobal *partition);
};

}



#endif /* SRC_WRAPPERS_WPARMETIS_H_ */
