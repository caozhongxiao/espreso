
#include "../../config/ecf/environment.h"

#include "../../physics/instance.h"

#include "../../wrappers/hypre/hypre.h"
#include "../../solver/generic/SparseMatrix.h"


#include "../../basis/utilities/utils.h"
#include "hypresolver.h"

using namespace espreso;

HYPRESolver::HYPRESolver(Instance *instance, const HypreConfiguration &configuration)
: _instance(instance), _hypreData(NULL)
{
	_configuration = configuration;
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
		_hypreData = new HypreData(environment->MPICommunicator, _instance->K[0].rows - _instance->K[0].haloRows);
	}

	if (matrices & Matrices::K) {
		size_t prefix = _instance->K[0].CSR_I_row_indices[_instance->K[0].haloRows] - 1;
		_hypreData->insertCSR(
				_instance->K[0].rows - _instance->K[0].haloRows, 0,
				_instance->K[0].CSR_I_row_indices.data() + _instance->K[0].haloRows,
				_instance->K[0].CSR_J_col_indices.data() + prefix,
				_instance->K[0].CSR_V_values.data() + prefix,
				_instance->f[0].data() + _instance->K[0].haloRows);
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
			_instance->K[0].rows - _instance->K[0].haloRows,
			_instance->primalSolution[0].data() + _instance->K[0].haloRows);
}

void HYPRESolver::finalize()
{

}
