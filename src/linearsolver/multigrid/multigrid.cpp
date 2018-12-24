
#include "multigrid.h"

#include "../../solver/generic/SparseMatrix.h"

#include "../../globals/run.h"
#include "../../basis/utilities/utils.h"
#include "../../physics/dataholder.h"
#include "../../config/ecf/solver/multigrid.h"
#include "../../wrappers/hypre/whypre.h"

using namespace espreso;

MultigridSolver::MultigridSolver(MultigridConfiguration &configuration)
: _configuration(configuration), _hypreData(NULL)
{

}

MultigridSolver::~MultigridSolver()
{
	if (_hypreData) {
		delete _hypreData;
	}
}

double& MultigridSolver::precision()
{
	return _configuration.precision;
}

void MultigridSolver::update(Matrices matrices)
{
	if (_hypreData == NULL) {
		_hypreData = new HypreData(run::data->K[0].rows - run::data->K[0].haloRows);
	}

	if (matrices & Matrices::K) {
		size_t prefix = run::data->K[0].CSR_I_row_indices[run::data->K[0].haloRows] - 1;
		_hypreData->insertCSR(
				run::data->K[0].rows - run::data->K[0].haloRows, 0,
				run::data->K[0].CSR_I_row_indices.data() + run::data->K[0].haloRows,
				run::data->K[0].CSR_J_col_indices.data() + prefix,
				run::data->K[0].CSR_V_values.data() + prefix,
				run::data->f[0].data() + run::data->K[0].haloRows);
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
void MultigridSolver::solve()
{
	HYPRE::solve(_configuration, *_hypreData,
			run::data->K[0].rows - run::data->K[0].haloRows,
			run::data->primalSolution[0].data() + run::data->K[0].haloRows);
}

void MultigridSolver::finalize()
{

}
