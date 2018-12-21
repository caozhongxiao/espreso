
#include "composer.h"

#include "../controllers/controller.h"

#include "../../dataholder.h"

#include "../../../globals/run.h"
#include "../../../basis/matrices/matrixtype.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../solver/generic/SparseMatrix.h"



using namespace espreso;

Composer::Composer(Controler &controler)
: _controler(controler), _DOFMap(NULL)
{

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
	run::storeSolution();
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
		std::fill(run::data->K[domain].CSR_V_values.begin(), run::data->K[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::M) {
		std::fill(run::data->M[domain].CSR_V_values.begin(), run::data->M[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::f) {
		std::fill(run::data->f[domain].begin(), run::data->f[domain].end(), 0);
	}
	if (matrices & Matrices::R) {
		std::fill(run::data->R[domain].begin(), run::data->R[domain].end(), 0);
	}
}

void Composer::keepK()
{
	#pragma omp parallel for
	for (size_t d = 0; d < run::mesh->elements->ndomains; d++) {
		run::data->origK[d] = run::data->K[d];
	}
}

void Composer::sum(NodeData *z, double alfa, NodeData* a, double beta, NodeData *b)
{

}

double Composer::multiply(NodeData *x, NodeData* y)
{
	return 0;
}
