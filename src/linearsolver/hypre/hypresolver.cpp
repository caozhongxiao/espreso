
#include "hypresolver.h"
#include "physics/assembler/dataholder.h"

#include "basis/utilities/utils.h"
#include "config/ecf/linearsolver/hypre/hypre.h"

#include "solver/generic/SparseMatrix.h"
#include "wrappers/hypre/hyprewrapper.h"

using namespace espreso;

HYPRESolver::HYPRESolver(DataHolder *data, HypreConfiguration &configuration)
: LinearSolver(data), _configuration(configuration), _hypreData(NULL)
{

}

HYPRESolver::~HYPRESolver()
{
	if (_hypreData) {
		delete _hypreData;
	}
}

double& HYPRESolver::precision()
{
	return _configuration.pcg.relative_conv_tol;
}

void HYPRESolver::update(Matrices matrices)
{
	if (_hypreData == NULL) {
		_hypreData = new HypreData(_data->K[0].rows - _data->K[0].haloRows);
	}

	if (matrices & Matrices::K) {
		size_t prefix = _data->K[0].CSR_I_row_indices[_data->K[0].haloRows] - 1;
		_hypreData->insertCSR(
				_data->K[0].rows - _data->K[0].haloRows, 0,
				_data->K[0].CSR_I_row_indices.data() + _data->K[0].haloRows,
				_data->K[0].CSR_J_col_indices.data() + prefix,
				_data->K[0].CSR_V_values.data() + prefix,
				_data->f[0].data() + _data->K[0].haloRows);
	}

//	if (matrices & Matrices::B1) {
//		_hypreData->insertIJV(
//				_instance->B1[0].rows, _instance->K[0].rows,
//				_instance->B1[0].nnz,
//				_instance->B1[0].I_row_indices.data(),
//				_instance->B1[0].J_col_indices.data(),
//				_instance->B1[0].V_values.data(),
//				_instance->B1c[0].data());
//	}

	_hypreData->finalizePattern();
}

// run solver and store primal and dual solution
void HYPRESolver::solve()
{
	HYPRE::solve(_configuration, *_hypreData,
			_data->K[0].rows - _data->K[0].haloRows,
			_data->primalSolution[0].data() + _data->K[0].haloRows);
}

void HYPRESolver::finalize()
{

}
