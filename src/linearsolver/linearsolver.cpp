
#include "physics/assembler/dataholder.h"
#include "linearsolver.h"

#include "basis/logging/logging.h"
#include "basis/utilities/utils.h"
#include "config/ecf/environment.h"

#include "solver/generic/SparseMatrix.h"

using namespace espreso;

template <typename TData>
static void store(TData &data, size_t domain, const std::string &name) {
	std::ofstream os(Logging::prepareFile(domain, name));
	os.precision(10);
	os << data;
	os.close();
};

LinearSolver::LinearSolver(DataHolder *data)
: _data(data)
{

}

void LinearSolver::solve(Matrices matrices)
{
	update(matrices);

	if (environment->print_matrices) {
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::BLUE << "STORE ASSEMBLED LINEAR SYSTEM";
		for (size_t d = 0; d < _data->K.size(); d++) {
			store(_data->K, d, "K");
			store(_data->N1, d, "N1");
			store(_data->N2, d, "N2");
			store(_data->RegMat, d, "RegMat");
			store(_data->M, d, "M");
			store(_data->R, d, "R");
			store(_data->f, d, "f");
			store(_data->B0, d, "B0");
			store(_data->B1, d, "B1");
			store(_data->B1c, d, "B1c");
			store(_data->B1duplicity, d, "B1duplicity");
			store(_data->primalSolution, d, "solution");
			store(_data->dualSolution, d, "dualSolution");
		}
	}

	solve();

	if (environment->print_matrices) {
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::BLUE << "STORE ASSEMBLED SYSTEM SOLUTION";
		for (size_t d = 0; d < _data->K.size(); d++) {
			store(_data->primalSolution, d, "solution");
			store(_data->dualSolution, d, "dualSolution");
		}
	}
}



