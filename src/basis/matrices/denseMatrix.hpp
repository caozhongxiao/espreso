#include "denseMatrix.h"

#include <map>

namespace espreso {

template<typename Tindices>
DenseMatrix::DenseMatrix(const SparseCSRMatrix<Tindices> &other): Matrix(other.rows(), other.columns(), DenseMatrixIndexing)
{
	std::vector<double>(rows() * columns(), 0).swap(_values);

	const Tindices *rowPtrs = other.rowPtrs();
	const Tindices *columnIndices = other.columnIndices();
	const double *values = other.values();

	for (Tindices r = 0; r < _rows; r++) {
		for (Tindices i = rowPtrs[r]; i < rowPtrs[r + 1]; i++) {
			set(r, columnIndices[i - other.indexing()] - other.indexing(), values[i - other.indexing()]);
		}
	}
}


template<typename Tindices>
DenseMatrix& DenseMatrix::operator=(const SparseCSRMatrix<Tindices> &other)
{
	DenseMatrix tmp(other);
	assign(*this, tmp);
	return *this;
}

}
