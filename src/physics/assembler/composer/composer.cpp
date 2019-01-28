
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "physics/assembler/dataholder.h"
#include "composer.h"

#include "physics/assembler/controllers/controller.h"

#include "basis/matrices/matrixtype.h"
#include "basis/utilities/communication.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "solver/generic/SparseMatrix.h"


using namespace espreso;

Composer::Composer(Controller &controler)
: _controler(controler), _foreignDOFs(0), _DOFMap(NULL)
{
	data = new DataHolder();
}

Composer::~Composer()
{
	delete data;
}

void Composer::initData()
{
	_controler.initData();
}

void Composer::nextTime()
{
	_controler.nextTime();
}

void Composer::parametersChanged()
{
	_controler.parametersChanged();
}

void Composer::processSolution()
{
	_controler.processSolution();
	info::mesh->storeSolution();
}

void Composer::insertKPattern(IJ *target, esint *begin, esint *end, MatrixType mtype)
{
	switch (mtype) {
	case MatrixType::REAL_SYMMETRIC_INDEFINITE:
	case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		for (auto row = begin, colbegin = begin; row != end; ++row, ++colbegin) {
			for (auto col = colbegin; col != end; ++col, ++target) {
				if (*row <= *col) {
					target->row = *row;
					target->column = *col;
				} else {
					target->row = *col;
					target->column = *row;
				}
			}
		}
		break;
	case MatrixType::REAL_UNSYMMETRIC:
		for (auto row = begin; row != end; ++row) {
			for (auto col = begin; col != end; ++col, ++target) {
				target->row = *row;
				target->column = *col;
			}
		}
		break;
	}
}

void Composer::clearMatrices(Matrices matrices, esint domain)
{
	if (matrices & Matrices::K) {
		std::fill(data->K[domain].CSR_V_values.begin(), data->K[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::M) {
		std::fill(data->M[domain].CSR_V_values.begin(), data->M[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::f) {
		std::fill(data->f[domain].begin(), data->f[domain].end(), 0);
	}
	if (matrices & Matrices::R) {
		std::fill(data->R[domain].begin(), data->R[domain].end(), 0);
	}
}

void Composer::keepK()
{
	data->origK.resize(data->K.size());
	#pragma omp parallel for
	for (size_t d = 0; d < data->K.size(); d++) {
		data->origK[d] = data->K[d];
	}
}

void Composer::keepRHS()
{
	data->origF.resize(data->f.size());
	#pragma omp parallel for
	for (size_t d = 0; d < data->f.size(); d++) {
		data->origF[d] = data->f[d];
	}
}

void Composer::keepSolverRHS()
{
	data->solverF.resize(data->f.size());
	#pragma omp parallel for
	for (size_t d = 0; d < data->f.size(); d++) {
		data->solverF[d] = data->f[d];
	}
}

void Composer::enrichSolution(double alfa, NodeData* x)
{
	for (size_t i = 0; i < _controler.solution()->data.size(); i++) {
		_controler.solution()->data[i] += alfa * x->data[i];
	}
}

void Composer::RHSMinusR()
{
	#pragma omp parallel for
	for (size_t d = 0; d < data->f.size(); d++) {
		for (size_t i = 0; i < data->f[d].size(); i++) {
			data->f[d][i] -= data->R[d][i];
		}
	}
}

void Composer::applyOriginalK(NodeData* result, NodeData* x)
{
	apply(data->origK, result->data, x->data);
}

void Composer::applyM(NodeData* result, NodeData* x)
{
	apply(data->M, result->data, x->data);
}

void Composer::sum(NodeData *z, double alfa, NodeData* a, double beta, NodeData *b)
{
	for (size_t i = 0; i < z->data.size(); i++) {
		z->data[i] = alfa * a->data[i] + beta * b->data[i];
	}
}

double Composer::norm(NodeData *x, NodeData* y)
{
	esint foreignPrefix = info::mesh->nodes->size - info::mesh->nodes->uniqueSize;
	double square = 0;
	for (size_t i = foreignPrefix * x->dimension; i < x->data.size(); ++i) {
		square += x->data[i] * y->data[i];
	}

	double sum = 0;
	MPI_Allreduce(&square, &sum, 1, MPI_DOUBLE, MPI_SUM, info::mpi::comm);

	return std::sqrt(sum);
}
