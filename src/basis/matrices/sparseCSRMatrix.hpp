#include "sparseCSRMatrix.h"
#include "mkl.h"

namespace espreso {

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(): Matrix(CSRMatrixIndexing)
{
	_rowPtrs.assign(2, _indexing);
}

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(size_t rows, size_t columns): Matrix(rows, columns, CSRMatrixIndexing)
{
	_rowPtrs.assign(rows + 1, _indexing);
}

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(const DenseMatrix &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	Tindices nnz = other.nonZeroValues();
	Tindices rows = _rows;
	Tindices columns = _columns;
	_rowPtrs.resize(other.rows() + 1);
	_columnIndices.resize(nnz);
	_values.resize(nnz);

	Tindices info;
	Tindices job[6] = {
		0,					// convert from dense to CSR
		other.indexing(),	// indexing of dense matrix
		indexing(),			// indexing of CSR matrix
		2,					// full matrix
		nnz,				// number of non-zero values
		1					// generate full output
	};

	mkl_ddnscsr (
		job, &rows, &columns,
		const_cast<double*>(other.values()), &columns,
		values(), columnIndices(), rowPtrs(),
		&info);
}

template<typename Tindices>
SparseCSRMatrix<Tindices>::SparseCSRMatrix(SparseVVPMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), CSRMatrixIndexing)
{
	other.shrink();
	size_t nnz = other.nonZeroValues();
	_rowPtrs.reserve(other.rows() + 1);
	_columnIndices.reserve(nnz);
	_values.reserve(nnz);

	const VVP<Tindices> &values = other.values();

	_rowPtrs.push_back(_indexing);
	for (size_t row = 0; row < _rows; row++) {
		for (size_t column = 0; column < values[row].size(); column++) {
			_columnIndices.push_back(values[row][column].first + _indexing);
			_values.push_back(values[row][column].second);
		}
		_rowPtrs.push_back(_values.size() + _indexing);
	}
}

template<typename Tindices>
SparseCSRMatrix<Tindices>& SparseCSRMatrix<Tindices>::operator=(const DenseMatrix &other)
{
	SparseCSRMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
SparseCSRMatrix<Tindices>& SparseCSRMatrix<Tindices>::operator=(SparseVVPMatrix<Tindices> &other)
{
	SparseCSRMatrix<Tindices> tmp(other);
	assign(*this, tmp);
	return *this;
}

template<typename Tindices>
void SparseCSRMatrix<Tindices>::multiply(SparseCSRMatrix<Tindices> &A, SparseCSRMatrix<Tindices> &B, bool transposeA)
{
	if (_indexing == Matrix::ZeroBased) {
		ESINFO(ERROR) << "Multiplication of two CSR matrices with zero based indexing is not supported.";
	}
	if (transposeA) {
		if (A.rows() != B.rows()) {
			ESINFO(ERROR) << "Matrix multiplication: matrices have incorrect dimensions.";
		}
	} else {
		if (A.columns() != B.rows()) {
			ESINFO(ERROR) << "Matrix multiplication: matrices have incorrect dimensions.";
		}
	}

	Tindices request = 0;
	Tindices sort = 8;	// C is sorted
	Tindices m = A.rows();
	Tindices n = A.columns();
	Tindices k = B.columns();
	Tindices nnz;		// not used
	Tindices info;

	_rows = transposeA ? A.columns() : A.rows();;
	_columns = k;
	_rowPtrs.resize(_rows + 1);

	request = 1;	// compute only rowPrts
	mkl_dcsrmultcsr(
		transposeA ? "t" : "n",
		&request,
		&sort,
		&m, &n, &k,
		A.values(), A.columnIndices(), A.rowPtrs(),
		B.values(), B.columnIndices(), B.rowPtrs(),
		values(), columnIndices(), rowPtrs(),
		&nnz,
		&info);

	_columnIndices.resize(_rowPtrs.back() - 1);
	_values.resize(_rowPtrs.back() - 1);

	request = 2;	// compute the rest of the matrix
	mkl_dcsrmultcsr(
		transposeA ? "t" : "n",
		&request,
		&sort,
		&m, &n, &k,
		A.values(), A.columnIndices(), A.rowPtrs(),
		B.values(), B.columnIndices(), B.rowPtrs(),
		values(), columnIndices(), rowPtrs(),
		&nnz,
		&info);
}

template<typename Tindices>
void SparseCSRMatrix<Tindices>::resize(size_t rows, size_t columns)
{
	_rows = rows;
	_columns = columns;
	_rowPtrs.resize(rows + 1, _rowPtrs.back());
}

template<typename Tindices>
void SparseCSRMatrix<Tindices>::transpose()
{
	Tindices job[6] = {
			0,		// CSR to CSC
			_indexing,
			_indexing,
			0,
			0,
			1
	};

	size_t size;
	if (_rows < _columns) {
		size = _columns;
		_rowPtrs.resize(size + 1, _rowPtrs.back());
	} else {
		size = _rows;
	}

	std::vector<Tindices> colPtrs(size + 1);
	std::vector<Tindices> rowIndices(_columnIndices.size());
	std::vector<double> vals(_values.size());

	Tindices n = size;
	Tindices info;

	mkl_dcsrcsc(
			job, &n,
			values(), columnIndices(), rowPtrs(),
			&vals[0], &rowIndices[0], &colPtrs[0],
			&info);

	colPtrs.resize(_columns + 1);
	_rowPtrs.swap(colPtrs);
	_columnIndices.swap(rowIndices);
	_values.swap(vals);
	size_t tmp = _rows;
	_rows = _columns;
	_columns = tmp;
}

template<typename Tindices>
int SparseCSRMatrix<Tindices>::gmresSolve(double *rhs, double *computed_solution, double tolerance, int maxiter) {
	//---------------------------------------------------------------------------
	// Define arrays for the coefficient matrix
	// Compressed sparse row storage is used for sparse representation
	//---------------------------------------------------------------------------
	//TODO:Check MKL_INT length and compare with Tindices
	int size = 128;
	Tindices *ia = &_rowPtrs[0];
	Tindices *ja = &_columnIndices[0];
	double *A = &_values[0];

	int N = columns();

	//---------------------------------------------------------------------------
	// Allocate storage for the ?par parameters and the solution/rhs/residual vectors
	//---------------------------------------------------------------------------
	Tindices ipar[size];
	double dpar[size];
	std::vector<double> tmp(N * (2 * N + 1)+(N * (N + 9)) / 2 + 1);

	//---------------------------------------------------------------------------
	// Some additional variables to use with the RCI (P)FGMRES solver
	//---------------------------------------------------------------------------
	Tindices itercount, ierr = 0;
	Tindices RCI_request, ivar;
	double dvar;
	char cvar, cvar1, cvar2;

	ivar = N;
	cvar = 'N';
	//---------------------------------------------------------------------------
	// Initialize the initial guess
	//---------------------------------------------------------------------------
	for (int i = 0; i < N; i++) {
		computed_solution[i] = 0.0;
	}

	//---------------------------------------------------------------------------
	// Initialize the solver
	//---------------------------------------------------------------------------
	dfgmres_init(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, &tmp[0]);
	if (RCI_request != 0) return -1;

	//---------------------------------------------------------------------------
	// Set the desired parameters:
	// https://software.intel.com/en-us/node/521710
	//---------------------------------------------------------------------------
	ipar[4] = maxiter;
	dpar[0] = tolerance;

	ipar[7] = 1;
	ipar[8] = 1;
	ipar[9] = 0;
	ipar[10] = 0;
	ipar[11] = 1;

	//---------------------------------------------------------------------------
	// Check the correctness and consistency of the newly set parameters
	//---------------------------------------------------------------------------
	dfgmres_check(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, &tmp[0]);
	if (RCI_request != 0) return -1;

	//---------------------------------------------------------------------------
	// Compute the solution by RCI (P)FGMRES solver with preconditioning
	// Reverse Communication starts here
	//---------------------------------------------------------------------------
	while (true) {
		dfgmres(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, &tmp[0]);
		//---------------------------------------------------------------------------
		// If RCI_request=0, then the solution was found with the required precision
		//---------------------------------------------------------------------------
		//std::cout<<"RCI "<<RCI_request<<std::endl;

		if (RCI_request == 0) break;
		//---------------------------------------------------------------------------
		// If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
		// and put the result in vector tmp[ipar[22]-1]
		//---------------------------------------------------------------------------
		if (RCI_request == 1) {
			mkl_dcsrgemv(&cvar, &ivar, A, ia, ja, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);
			continue;
		}

		break;
	}
	//---------------------------------------------------------------------------
	// Reverse Communication ends here
	// Get the current iteration number and the FGMRES solution (DO NOT FORGET to
	// call dfgmres_get routine as computed_solution is still containing
	// the initial guess!). Request to dfgmres_get to put the solution
	// into vector computed_solution[N] via ipar[12]
	//---------------------------------------------------------------------------
	ipar[12] = 0;
	dfgmres_get(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, &tmp[0], &itercount);
	//---------------------------------------------------------------------------
	// Print solution vector: computed_solution[N] and the number of iterations: itercount
	//---------------------------------------------------------------------------
	//printf("The system has been solved \n");
	//printf("\nNumber of iterations: %d\n" , itercount);
	//for (int i = 0; i < N; i++) {
	//	printf("%f ",computed_solution[i]);
	//}
	//printf("\n");
	return 0;
}

}

