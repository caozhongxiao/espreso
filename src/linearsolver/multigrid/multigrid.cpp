
#include "../../basis/utilities/utils.h"
#include "../../wrappers/hypre/hypre.h"
#include "../../solver/generic/SparseMatrix.h"
#include "../../config/ecf/environment.h"
#include "multigrid.h"

using namespace espreso;

MultigridSolver::MultigridSolver(Instance *instance, const MultigridConfiguration &configuration)
: instance(instance),
  timeEvalMain("ESPRESO Multigrid Solver Overal Timing")
{
	this->configuration = configuration;
}

MultigridSolver::~MultigridSolver() {
	// TODO Auto-generated destructor stub
}


// make partial initialization according to updated matrices
void MultigridSolver::update(Matrices matrices)
{

}


// run solver and store primal and dual solution
void MultigridSolver::solve()
{
	SparseMatrix& k = instance->K[0];

	HypreRegion region;
	region.nrows = k.rows;
	for (eslocal i = 0; i < k.rows; i++) {
		region.ncols.push_back(k.CSR_I_row_indices[i + 1]-k.CSR_I_row_indices[i]);
		region.rows.push_back(i+1);
	}
	region.cols = k.CSR_J_col_indices.data();
	region.values = k.CSR_V_values.data();

	std::vector<HypreRegion> values;
	values.push_back(region);

	//TODO: repair this
	instance->f[0][0]=1;
	instance->f[0][15]=2;

	HYPRE::Solve(this->configuration,
			environment->MPIrank, environment->MPIsize, environment->MPICommunicator,
			1, instance->K[0].rows, values,
			instance->f[0].data(),
			instance->primalSolution[0].data());

}

void MultigridSolver::finalize()
{
}
