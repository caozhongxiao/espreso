
#ifndef SRC_WRAPPERS_METIS_WMETIS_H_
#define SRC_WRAPPERS_METIS_WMETIS_H_

namespace espreso {

struct METISConfiguration;

struct METIS {

	static esint call(
			const METISConfiguration &options,
			esint verticesCount,
			esint *eframes, esint *eneighbors,
			esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
			esint parts, esint *partition);

};

}



#endif /* SRC_WRAPPERS_METIS_WMETIS_H_ */
