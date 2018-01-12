
#ifndef SRC_WRAPPERS_MATH_MATH_H_
#define SRC_WRAPPERS_MATH_MATH_H_

namespace espreso {

struct MATH {

	static void upCSRMatVecProduct(eslocal rows, eslocal cols, eslocal *mRows, eslocal *mCols, float *mVals, float *vVals, float *result);
	static void upCSRMatVecProduct(eslocal rows, eslocal cols, eslocal *mRows, eslocal *mCols, double *mVals, double *vVals, double *result);

	static void vecScale(eslocal size, float alpha, float *vVals);
	static void vecScale(eslocal size, double alpha, double *vVals);

	static double vecNorm(eslocal size, float *vVals);
	static double vecNorm(eslocal size, double *vVals);

	static eslocal vecNormMaxIndex(eslocal size, float *vVals);
	static eslocal vecNormMaxIndex(eslocal size, double *vVals);
};

}



#endif /* SRC_WRAPPERS_MATH_MATH_H_ */
