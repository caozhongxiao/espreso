
#ifndef SRC_WRAPPERS_METIS_WMETIS_H_
#define SRC_WRAPPERS_METIS_WMETIS_H_

namespace espreso {

struct METIS {

	static eslocal call(
			eslocal verticesCount,
			eslocal *eframes, eslocal *eneighbors,
			eslocal verticesWeightCount, eslocal *verticesWeights, eslocal *edgeWeights,
			eslocal parts, eslocal *partition);

};

}



#endif /* SRC_WRAPPERS_METIS_WMETIS_H_ */
