
#ifndef SRC_WRAPPERS_WMETIS_H_
#define SRC_WRAPPERS_WMETIS_H_

namespace espreso {

struct METIS {

	static esglobal call(
			esglobal verticesCount,
			esglobal *eframes, esglobal *eneighbors,
			esglobal verticesWeightCount, esglobal *verticesWeights, esglobal *edgeWeights,
			esglobal parts, esglobal *partition);

};

}



#endif /* SRC_WRAPPERS_WMETIS_H_ */
