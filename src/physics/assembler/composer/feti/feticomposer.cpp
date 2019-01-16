
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
