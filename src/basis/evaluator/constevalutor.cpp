
#include "constevaluator.h"

using namespace espreso;

void ConstEvaluator::evaluate(eslocal size, eslocal increment, const Point* cbegin, const double* tbegin, double time, double *results) const
{
	for (eslocal i = 0; i < size; ++i) {
		results[i * increment] = _value;
	}
}


