
#ifndef SRC_WRAPPERS_MATH_MATH_H_
#define SRC_WRAPPERS_MATH_MATH_H_

namespace espreso {

struct MATH {

	static void setNumberOfThreads(int numberOfThreads);

	static void upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, float *mVals, float *vVals, float *result);
	static void upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, double *mVals, double *vVals, double *result);

	static void vecScale(esint size, float alpha, float *vVals);
	static void vecScale(esint size, double alpha, double *vVals);

	static double vecNorm(esint size, float *vVals);
	static double vecNorm(esint size, double *vVals);

	static esint vecNormMaxIndex(esint size, float *vVals);
	static esint vecNormMaxIndex(esint size, double *vVals);

	// C = alpha * A * B + beta * C
	static void DenseMatDenseMatRowMajorProduct(
			double alpha, bool transposeA, esint aRows, esint aCols, double* aVals,
			bool transposeB, esint bCols, double* bVals,
			double beta, double* cVals);

	struct SOLVER {

		static void GMRESUpCRSMat(
				esint rows, esint cols, esint *mRows, esint *mCols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, esint maxIterations, esint &itercount);

		static void GMRESDenseRowMajorMat(
				esint rows, esint cols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, esint maxIterations, esint &itercount);

		static void GMRESUpperSymetricColumnMajorMat(
				esint cols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, esint maxIterations, esint &itercount);

		static esint directUpperSymetricIndefiniteColumnMajor(
				esint cols, double *m_packed_values,
				esint nrhs, double *rhsVals);
	};
};

}



#endif /* SRC_WRAPPERS_MATH_MATH_H_ */
