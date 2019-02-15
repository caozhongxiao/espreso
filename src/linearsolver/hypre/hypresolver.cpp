
#include "hypresolver.h"
#include "physics/assembler/dataholder.h"

#include "basis/logging/logging.h"
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
	switch (_configuration.solver_type) {
	case HypreConfiguration::SOLVER_TYPE::BiCGSTAB:
		return _configuration.bicgstab.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::BoomerAMG:
		return _configuration.boomeramg.convergence_tolerance;
	case HypreConfiguration::SOLVER_TYPE::CGNR:
		return _configuration.cgnr.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::FlexGMRES:
		return _configuration.flexgmres.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::GMRES:
		return _configuration.gmres.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::LGMRES:
		return _configuration.lgmres.relative_conv_tol;
	case HypreConfiguration::SOLVER_TYPE::PCG:
		return _configuration.pcg.relative_conv_tol;
	default:
		ESINFO(GLOBAL_ERROR) << "Required precision of unknown solver.";
		exit(0);
	}
}

void HYPRESolver::update(Matrices matrices)
{
	if (_hypreData == NULL) {
		_hypreData = new HypreData(_data->K[0].rows - _data->K[0].haloRows);
	}

	if (matrices & (Matrices::K | Matrices::f)) {
		size_t prefix = _data->K[0].CSR_I_row_indices[_data->K[0].haloRows] - 1;
		_hypreData->insertCSR(
				_data->K[0].rows - _data->K[0].haloRows, 0,
				_data->K[0].CSR_I_row_indices.data() + _data->K[0].haloRows,
				_data->K[0].CSR_J_col_indices.data() + prefix,
				_data->K[0].CSR_V_values.data() + prefix,
				_data->f[0].data() + _data->K[0].haloRows);
		_hypreData->finalizePattern();
	}
}

// run solver and store primal and dual solution
void HYPRESolver::solve()
{
	_hypreData->solve(_configuration,
			_data->K[0].rows - _data->K[0].haloRows,
			_data->primalSolution[0].data() + _data->K[0].haloRows);
}

void HYPRESolver::finalize()
{

}
