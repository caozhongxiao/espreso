
#include "wrappers/math/math.h"
#include "esinfo/mpiinfo.h"

#include "mkl.h"

using namespace espreso;

void MATH::setNumberOfThreads(int numberOfThreads) {
	mkl_set_num_threads(numberOfThreads);
}

void MATH::upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, float *mVals, float *vVals, float *result)
{
	mkl_cspblas_scsrsymv("U", &rows, mVals, mRows, mCols, vVals, result);
}

void MATH::upCSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, double *mVals, double *vVals, double *result)
{
	mkl_cspblas_dcsrsymv("U", &rows, mVals, mRows, mCols, vVals, result);
}

void MATH::vecScale(esint size, float alpha, float *vVals)
{
	esint incr = 1;
	cblas_sscal(size, alpha, vVals, incr);
}

void MATH::vecScale(esint size, double alpha, double *vVals)
{
	esint incr = 1;
	cblas_dscal(size, alpha, vVals, incr);
}

double MATH::vecNorm(esint size, float *vVals)
{
	esint incr = 1;
	return snrm2(&size, vVals, &incr);
}

double MATH::vecNorm(esint size, double *vVals)
{
	esint incr = 1;
	return dnrm2(&size, vVals, &incr);
}

esint MATH::vecNormMaxIndex(esint size, float *vVals)
{
	esint incr = 1;
	return cblas_isamax(size, vVals, incr);
}

esint MATH::vecNormMaxIndex(esint size, double *vVals)
{
	esint incr = 1;
	return cblas_idamax(size, vVals, incr);
}

void MATH::DenseMatDenseMatRowMajorProduct(
			double alpha, bool transposeA, esint aRows, esint aCols, double* aVals,
			bool transposeB, esint bCols, double* bVals,
			double beta, double* cVals)
{
	cblas_dgemm(CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans, transposeB ? CblasTrans : CblasNoTrans,
			aRows, bCols, aCols, alpha, aVals, aRows, bVals, aCols, beta, cVals, aRows);
}

void MATH::CSRMatVecProduct(esint rows, esint cols, esint *mRows, esint *mCols, double *mVals, double *vVals, double *result)
{
	mkl_dcsrgemv("N", &rows, mVals, mRows, mCols, vVals, result);
}
