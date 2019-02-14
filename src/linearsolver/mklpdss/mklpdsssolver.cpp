
#include "mklpdsssolver.h"

#include "config/ecf/linearsolver/mklpdss.h"
#include "physics/assembler/dataholder.h"

#include "solver/generic/SparseMatrix.h"
#include "wrappers/mklpdss/mklpdsswrapper.h"

using namespace espreso;

MKLPDSSSolver::MKLPDSSSolver(DataHolder *data, MKLPDSSConfiguration &configuration)
: LinearSolver(data), _configuration(configuration), _mklpdssData(NULL), _precision(1e-12)
{

}

MKLPDSSSolver::~MKLPDSSSolver()
{
	if (_mklpdssData) {
		delete _mklpdssData;
	}
}

void MKLPDSSSolver::update(Matrices matrices)
{
	if (_mklpdssData == NULL) {
		_mklpdssData = new MKLPDSSData(_data->K[0].rows - _data->K[0].haloRows);
	}
	if (matrices & Matrices::K) {
		size_t prefix = _data->K[0].CSR_I_row_indices[_data->K[0].haloRows] - 1;
		_mklpdssData->insertCSR(
				_data->K[0].mtype,
				_data->K[0].CSR_I_row_indices.data() + _data->K[0].haloRows,
				_data->K[0].CSR_J_col_indices.data() + prefix,
				_data->K[0].CSR_V_values.data() + prefix,
				_data->f[0].data() + _data->K[0].haloRows);
	}
}

void MKLPDSSSolver::solve()
{
	_mklpdssData->solve(_configuration, _data->primalSolution[0].data() + _data->K[0].haloRows);
}

void MKLPDSSSolver::finalize()
{

}
