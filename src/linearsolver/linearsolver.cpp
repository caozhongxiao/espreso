
#include "linearsolver.h"

#include "globals/run.h"
#include "basis/logging/logging.h"
#include "basis/utilities/utils.h"
#include "config/ecf/environment.h"

#include "physics/dataholder.h"
#include "solver/generic/SparseMatrix.h"

using namespace espreso;

template <typename TData>
static void store(TData &data, size_t domain, const std::string &name) {
	std::ofstream os(Logging::prepareFile(domain, name));
	os.precision(10);
	os << data;
	os.close();
};

void LinearSolver::solve(Matrices matrices)
{
	update(matrices);

	if (environment->print_matrices) {
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::BLUE << "STORE ASSEMBLED LINEAR SYSTEM";
		for (size_t d = 0; d < run::data->K.size(); d++) {
			store(run::data->K, d, "K");
			store(run::data->N1, d, "N1");
			store(run::data->N2, d, "N2");
			store(run::data->RegMat, d, "RegMat");
			store(run::data->M, d, "M");
			store(run::data->R, d, "R");
			store(run::data->f, d, "f");
			store(run::data->B0, d, "B0");
			store(run::data->B1, d, "B1");
			store(run::data->B1c, d, "B1c");
			store(run::data->B1duplicity, d, "B1duplicity");
			store(run::data->primalSolution, d, "solution");
			store(run::data->dualSolution, d, "dualSolution");
		}
	}

	solve();

	if (environment->print_matrices) {
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::BLUE << "STORE ASSEMBLED SYSTEM SOLUTION";
		for (size_t d = 0; d < run::data->K.size(); d++) {
			store(run::data->primalSolution, d, "solution");
			store(run::data->dualSolution, d, "dualSolution");
		}
	}
}



