
#include "physics/assembler/dataholder.h"
#include "globalcomposer.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "solver/generic/SparseMatrix.h"

using namespace espreso;

NodeData* GlobalComposer::RHS()
{
	return NULL;
}

void GlobalComposer::KplusAlfaM(double alfa)
{
	data->K[0].MatAddInPlace(data->M[0], 'N', alfa);
}

void GlobalComposer::apply(std::vector<SparseMatrix> &matrices, NodeData *result, NodeData *x)
{
	matrices[0].MatVec(x->data, result->data, 'N', 0, 0, 0);
}

void GlobalComposer::enrichRHS(double alfa, NodeData* a)
{
	for (size_t i = 0; i < data->f[0].size(); ++i) {
		data->f[0][i] += alfa * a->data[i];
	}
}

void GlobalComposer::RHSMinusR()
{

}

void GlobalComposer::DirichletMinusRHS()
{

}

double GlobalComposer::residualNorm()
{
	return std::sqrt(0);
}

