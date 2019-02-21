
#include "mklpdsswrapper.h"

#include "esinfo/mpiinfo.h"

#include "basis/matrices/matrixtype.h"
#include "basis/utilities/communication.h"
#include "esinfo/eslog.h"

#include <vector>

#ifdef HAVE_MKL
#include "mkl_cluster_sparse_solver.h"

namespace espreso {
struct MKLPDSSDataHolder {
	void *pt[64];

	esint maxfct;
	esint mnum;

	esint mtype;
	esint phase;
	esint perm;
	esint n;
	esint nrhs;
	esint iparm[64];
	esint msglvl;
	int   comm;
	esint error;

	std::vector<esint> rowData, colData;
	std::vector<double> valData;

	esint *rowPtrs;
	esint *colIndices;
	double *values;
	double *rhsValues;
	double *solution;
};
}
#endif

using namespace espreso;

MKLPDSSData::MKLPDSSData(esint nrows)
: _roffset(nrows), _nrows(nrows), _data(NULL)
{
#ifdef HAVE_MKL
	_data = new MKLPDSSDataHolder();
	_data->n = Communication::exscan(_roffset);

	_data->maxfct = 1; // dummy
	_data->mnum = 1; // dummy

	_data->nrhs = 1;

	std::fill(_data->iparm, _data->iparm + 64, 0);
	_data->iparm[0] = 1; // Use filled values.
	// Fill-in reducing ordering for the input matrix.
	_data->iparm[1] = info::mpi::size > 1 ? 10 : 3; // MPI or parallel
	// Matrix input format.
	_data->iparm[39] = 2; // distributed A, x, rhs
	_data->iparm[40] = _roffset + 1;
	_data->iparm[41] = _roffset + _nrows;

	_data->msglvl = 0;
	_data->comm = MPI_Comm_c2f(info::mpi::comm);
#endif
}

MKLPDSSData::~MKLPDSSData()
{
#ifdef HAVE_MKL
	call(-1);
	delete _data;
#endif
}

void MKLPDSSData::insertK(MatrixType mtype, esint *rowPtrs, esint *colIndices, double *values)
{
#ifdef HAVE_MKL
	switch (mtype) {
	case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		_data->mtype = 2; break;
	case MatrixType::REAL_SYMMETRIC_INDEFINITE:
		_data->mtype = -2; break;
	case MatrixType::REAL_UNSYMMETRIC:
		_data->mtype = 1; break;
	}

	const bool fillStructure = !_data->colData.size();

	// pick only upper triangle (since composer does not set correct dirichlet in symmetric matrices)
	if (mtype == MatrixType::REAL_SYMMETRIC_INDEFINITE || mtype == MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		if (fillStructure) {
			_data->rowData.clear();
			_data->colData.clear();
			_data->rowData.reserve(_nrows + 1);
			_data->rowData.push_back(1);
		}
		_data->valData.clear();
		for (esint i = 0; i < _nrows; i++) {
			for (esint c = rowPtrs[i] - rowPtrs[0]; c < rowPtrs[i + 1] - rowPtrs[0]; c++) {
				if (_roffset + i <= colIndices[c] - 1) {
					if (fillStructure) {
						_data->colData.push_back(colIndices[c]);
					}
					_data->valData.push_back(values[c]);
				}
			}
			if (fillStructure) {
				_data->rowData.push_back(_data->colData.size() + 1);
			}
		}
		_data->colIndices = _data->colData.data();
		_data->values = _data->valData.data();
	} else {
		// row data have to be always renumbered
		if (!_data->rowData.size()) {
			_data->rowData.resize(_nrows + 1);
			for (size_t i = 0; i < _data->rowData.size(); i++) {
				_data->rowData[i] = rowPtrs[i] - rowPtrs[0] + 1;
			}
		}
		_data->colIndices = colIndices;
		_data->values = values;
	}

	_data->rowPtrs = _data->rowData.data();

	if (fillStructure) {
		call(12);
	} else {
		call(22);
	}
#endif
}

void MKLPDSSData::insertRHS(double *rhsValues)
{
#ifdef HAVE_MKL
	_data->rhsValues = rhsValues;
#endif
}

void MKLPDSSData::solve(const MKLPDSSConfiguration &configuration, double *solution)
{
#ifdef HAVE_MKL
	_data->solution = solution;
	call(33); // solve at once
	eslog::solver("MKL PDSS: Solved");
#endif
}

void MKLPDSSData::call(esint phase)
{
#ifdef HAVE_MKL
	_data->phase = phase;
	cluster_sparse_solver(
			_data->pt, &_data->maxfct, &_data->mnum,
			&_data->mtype,
			&_data->phase,
			&_data->n, _data->values, _data->rowPtrs, _data->colIndices,
			&_data->perm, &_data->nrhs, _data->iparm, &_data->msglvl,
			_data->rhsValues, _data->solution,
			&_data->comm, &_data->error);

	switch (_data->error) {
	case   0: break;
	case  -1: eslog::error("MKL PDSS: input inconsistent."); break;
	case  -2: eslog::error("MKL PDSS: not enough memory."); break;
	case  -3: eslog::error("MKL PDSS: reordering problem."); break;
	case  -4: eslog::error("MKL PDSS: zero pivot, numerical factorization or iterative refinement problem."); break;
	case  -5: eslog::error("MKL PDSS: unclassified (internal) error."); break;
	case  -6: eslog::error("MKL PDSS: reordering failed."); break;
	case  -7: eslog::error("MKL PDSS: diagonal matrix is singular."); break;
	case  -8: eslog::error("MKL PDSS: 32-bit integer overflow problem."); break;
	case  -9: eslog::error("MKL PDSS: not enough memory for OOC."); break;
	case -10: eslog::error("MKL PDSS: error opening OOC files."); break;
	case -11: eslog::error("MKL PDSS: read/write error with OOC files."); break;
	}

#endif
}

