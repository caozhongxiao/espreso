
#include "esinfo/meshinfo.h"
#include "physics/assembler/dataholder.h"
#include "composer.h"

#include "physics/assembler/controllers/controller.h"

#include "basis/matrices/matrixtype.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "solver/generic/SparseMatrix.h"


using namespace espreso;

Composer::Composer(Controler &controler)
: _controler(controler), _DOFMap(NULL)
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
	info::storeSolution();
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
	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		data->origK[d] = data->K[d];
	}
}

void Composer::sum(NodeData *z, double alfa, NodeData* a, double beta, NodeData *b)
{

}

double Composer::multiply(NodeData *x, NodeData* y)
{
	return 0;
}
