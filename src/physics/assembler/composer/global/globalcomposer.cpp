
#include "physics/assembler/dataholder.h"
#include "globalcomposer.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
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

void GlobalComposer::applyM(NodeData *y, NodeData *x)
{

}

void GlobalComposer::applyOriginalK(NodeData *y, NodeData *x)
{

}

void GlobalComposer::enrichRHS(double alfa, NodeData* a)
{

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

