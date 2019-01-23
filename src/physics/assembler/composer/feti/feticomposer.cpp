
#include "esinfo/meshinfo.h"
#include "physics/assembler/dataholder.h"
#include "feticomposer.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "solver/generic/SparseMatrix.h"

using namespace espreso;

NodeData* FETIComposer::RHS()
{
	return NULL;
}

void FETIComposer::KplusAlfaM(double alfa)
{
	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		data->K[d].MatAddInPlace(data->M[d], 'N', alfa);
	}
}

void FETIComposer::apply(std::vector<SparseMatrix> &matrices, NodeData *result, NodeData *x)
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

void FETIComposer::enrichRHS(double alfa, NodeData* x)
{
	std::vector<std::vector<double> > _x(data->f.size());

	for (size_t i = 0; i < data->f.size(); i++) {
		_x[i].resize(data->f[i].size());
	}
	divide(x, _x);

	#pragma omp parallel for
	for (size_t d = 0; d < data->f.size(); d++) {
		for (size_t i = 0; i < data->f[d].size(); ++i) {
			data->f[d][i] += _x[d][i];
		}
	}
}

void FETIComposer::RHSMinusR()
{

}

void FETIComposer::DirichletMinusRHS()
{

}

double FETIComposer::residualNorm()
{
	return std::sqrt(0);
}
