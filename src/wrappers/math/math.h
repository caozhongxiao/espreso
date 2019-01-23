
#ifndef SRC_WRAPPERS_MATH_MATH_H_
#define SRC_WRAPPERS_MATH_MATH_H_

namespace espreso {

namespace MATH {

	void setNumberOfThreads(int numberOfThreads);

	void upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, float *mVals, float *vVals, float *result);
	void upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, double *mVals, double *vVals, double *result);

	void vecScale(esint size, float alpha, float *vVals);
	void vecScale(esint size, double alpha, double *vVals);

	double vecNorm(esint size, float *vVals);
	double vecNorm(esint size, double *vVals);

	esint vecNormMaxIndex(esint size, float *vVals);
	esint vecNormMaxIndex(esint size, double *vVals);

	// C = alpha * A * B + beta * C
	void DenseMatDenseMatRowMajorProduct(
			double alpha, bool transposeA, esint aRows, esint aCols, double* aVals,
			bool transposeB, esint bCols, double* bVals,
			double beta, double* cVals);

	void CSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, double *mVals, double *vVals, double *result);

	namespace SOLVER {

		void GMRESUpCRSMat(
				esint rows, esint cols, esint *mRows, esint *mCols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, esint maxIterations, esint &itercount);

		void GMRESDenseRowMajorMat(
				esint rows, esint cols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, esint maxIterations, esint &itercount);

		void GMRESUpperSymetricColumnMajorMat(
				esint cols, double *mVals,
				double *rhsVals, double *results,
				double tolerance, esint maxIterations, esint &itercount);

		esint directUpperSymetricIndefiniteColumnMajor(
				esint cols, double *m_packed_values,
				esint nrhs, double *rhsVals);
	};
};

}



#endif /* SRC_WRAPPERS_MATH_MATH_H_ */
