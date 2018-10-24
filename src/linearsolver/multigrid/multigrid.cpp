
#include "multigrid.h"

#include "../../config/ecf/environment.h"

#include "../../assembler/instance.h"

#include "../../wrappers/hypre/hypre.h"
#include "../../solver/generic/SparseMatrix.h"


#include "../../basis/utilities/utils.h"

using namespace espreso;

MultigridSolver::MultigridSolver(Instance *instance, const MultigridConfiguration &configuration)
: _instance(instance), _hypreData(NULL)
{
	_configuration = configuration;
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
		_hypreData = new HypreData(environment->MPICommunicator, _instance->K[0].rows);
	}

	if (matrices & Matrices::K) {
		_hypreData->insertCSR(
				_instance->K[0].rows, 0,
				_instance->K[0].CSR_I_row_indices.data(),
				_instance->K[0].CSR_J_col_indices.data(),
				_instance->K[0].CSR_V_values.data(),
				_instance->f[0].data());
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
	HYPRE::solve(_configuration, *_hypreData, _instance->primalSolution[0].size(), _instance->primalSolution[0].data());
}

void MultigridSolver::finalize()
{

}
