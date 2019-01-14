
#include "feticomposer.h"

#include "physics/dataholder.h"

#include "globals/run.h"
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
	for (size_t d = 0; d < run::mesh->elements->ndomains; d++) {
		run::data->K[d].MatAddInPlace(run::data->M[d], 'N', alfa);
	}
}

void FETIComposer::applyM(NodeData *y, NodeData *x)
{

}

void FETIComposer::applyOriginalK(NodeData *y, NodeData *x)
{

}

void FETIComposer::enrichRHS(double alfa, NodeData* a)
{

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
