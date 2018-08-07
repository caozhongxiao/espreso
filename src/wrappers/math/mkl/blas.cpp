
#include "../math.h"

#include "mkl.h"

using namespace espreso;

void MATH::setNumberOfThreads(int numberOfThreads) {
	mkl_set_num_threads(numberOfThreads);
}

void MATH::upCSRMatVecProduct(eslocal rows, eslocal cols, eslocal *mRows, eslocal *mCols, float *mVals, float *vVals, float *result)
{
	mkl_cspblas_scsrsymv("U", &rows, mVals, mRows, mCols, vVals, result);
}

void MATH::upCSRMatVecProduct(eslocal rows, eslocal cols, eslocal *mRows, eslocal *mCols, double *mVals, double *vVals, double *result)
{
	mkl_cspblas_dcsrsymv("U", &rows, mVals, mRows, mCols, vVals, result);
}

void MATH::vecScale(eslocal size, float alpha, float *vVals)
{
	eslocal incr = 1;
	cblas_sscal(size, alpha, vVals, incr);
}

void MATH::vecScale(eslocal size, double alpha, double *vVals)
{
	eslocal incr = 1;
	cblas_dscal(size, alpha, vVals, incr);
}

double MATH::vecNorm(eslocal size, float *vVals)
{
	eslocal incr = 1;
	return snrm2(&size, vVals, &incr);
}

double MATH::vecNorm(eslocal size, double *vVals)
{
	eslocal incr = 1;
	return dnrm2(&size, vVals, &incr);
}

eslocal MATH::vecNormMaxIndex(eslocal size, float *vVals)
{
	eslocal incr = 1;
	return cblas_isamax(size, vVals, incr);
}

eslocal MATH::vecNormMaxIndex(eslocal size, double *vVals)
{
	eslocal incr = 1;
	return cblas_idamax(size, vVals, incr);
}

void MATH::DenseMatDenseMatRowMajorProduct(
			double alpha, bool transposeA, eslocal aRows, eslocal aCols, double* aVals,
			bool transposeB, eslocal bCols, double* bVals,
			double beta, double* cVals)
{
	cblas_dgemm(CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans, transposeB ? CblasTrans : CblasNoTrans,
			aRows, bCols, aCols, alpha, aVals, aRows, bVals, aCols, beta, cVals, aRows);
}
