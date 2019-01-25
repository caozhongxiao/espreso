
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "physics/assembler/dataholder.h"
#include "physics/assembler/controllers/controller.h"
#include "feticomposer.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "solver/generic/SparseMatrix.h"

using namespace espreso;

void FETIComposer::KplusAlfaM(double alfa)
{
	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		data->K[d].MatAddInPlace(data->M[d], 'N', alfa);
	}
}

void FETIComposer::apply(std::vector<SparseMatrix> &matrices, std::vector<double> &result, std::vector<double> &x)
{
	std::vector<std::vector<double> > _res(matrices.size()), _x(matrices.size());

	for (size_t i = 0; i < matrices.size(); i++) {
		_res[i].resize(matrices[i].rows);
		_x[i].resize(matrices[i].rows);
	}

	duply(x, _x);

	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		matrices[d].MatVec(_x[d], _res[d], 'N', 0, 0, 0);
	}

	gather(result, _res);
}

//void FETIComposer::DirichletMinusSolution()
//{
//	std::vector<std::vector<double> > solution;
//	duply(_controler.solution()->data, solution);
//	#pragma omp parallel for
//	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
//		const std::vector<esint> &ROWS = data->B1[d].I_row_indices;
//		for (size_t j = 0; j < ROWS.size() && ROWS[j] <= data->block[DataHolder::CONSTRAINT::DIRICHLET]; j++) {
//			data->B1c[d][j] -= solution[d][data->B1[d].J_col_indices[j] - 1];
//		}
//	}
//}

double FETIComposer::residualNormNumerator()
{
	double square = 0;
	std::vector<double> f, btlambda;
	gather(f, data->f);
	gather(btlambda, data->dualSolution);
	for (size_t i = _foreignDOFs; i < f.size(); i++) {
		square += (f[i] - btlambda[i]) * (f[i] - btlambda[i]);
	}

	double sum = 0;
	MPI_Allreduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, info::mpi::comm);

	return std::sqrt(sum);
}

double FETIComposer::residualNormDenominator()
{
	double square = 0;
	std::vector<double> f, r;
	gather(f, data->origF);
	gather(r, data->R);
	for (size_t i = _foreignDOFs, d = 0; i < f.size(); i++) {
		while (d < _dirichletMap.size() && (size_t)_dirichletMap[d] < i) { d++; }
		if (d == _dirichletMap.size() || (size_t)_dirichletMap[d] != i) {
			square += f[i] * f[i];
		} else {
			square += r[i] * r[i];
		}
	}

	double sum = 0;
	MPI_Allreduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, info::mpi::comm);

	return std::sqrt(sum);
}
