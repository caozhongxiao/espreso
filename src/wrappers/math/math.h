
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

	static void matMatProduct(
			eslocal aRows, eslocal aCols, double* aVals,
			eslocal bCols, double* bVals,
			double* cVals);

	struct SOLVER {

		static void GMRESUpCRSMat(
				eslocal rows, eslocal cols, eslocal *mRows, eslocal *mCols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, eslocal maxIterations);

		static void GMRESDenseRowMajorMat(
				eslocal rows, eslocal cols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, eslocal maxIterations);

		static void GMRESUpperSymetricColumnMajorMat(
				eslocal cols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, eslocal maxIterations);

		static int directUpperSymetricIndefiniteColumnMajor(
						eslocal cols, double *m_packed_values,
						eslocal nrhs, double *rhsVals);
	};
};

}



#endif /* SRC_WRAPPERS_MATH_MATH_H_ */
