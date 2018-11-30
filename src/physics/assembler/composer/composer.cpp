
#include "composer.h"

#include "../controllers/controller.h"

#include "../../instance.h"
#include "../../../basis/matrices/matrixtype.h"

#include "../../../solver/generic/SparseMatrix.h"


using namespace espreso;

Composer::Composer(Mesh &mesh, Step &step, Instance &instance, Controler &controler)
: _mesh(mesh), _step(step), _instance(instance), _controler(controler), _DOFMap(NULL)
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
}

void Composer::insertKPattern(IJ *target, eslocal *begin, eslocal *end, MatrixType mtype)
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

std::vector<double>& Composer::getSolutionStore()
{
	return _controler.getSolutionStore();
}

void Composer::clearMatrices(Matrices matrices, eslocal domain)
{
	if (matrices & Matrices::K) {
		std::fill(_instance.K[domain].CSR_V_values.begin(), _instance.K[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::M) {
		std::fill(_instance.M[domain].CSR_V_values.begin(), _instance.M[domain].CSR_V_values.end(), 0);
	}
	if (matrices & Matrices::f) {
		std::fill(_instance.f[domain].begin(), _instance.f[domain].end(), 0);
	}
	if (matrices & Matrices::R) {
		std::fill(_instance.R[domain].begin(), _instance.R[domain].end(), 0);
	}
}
