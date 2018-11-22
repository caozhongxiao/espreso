
#include "constevaluator.h"

using namespace espreso;

void ConstEvaluator::evalVector(eslocal size, eslocal increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	for (eslocal i = 0; i < size; ++i) {
		results[i * increment] = _value;
	}
}


void ConstEvaluator::evalFiltered(eslocal size, eslocal increment, eslocal *elements, eslocal *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	for (eslocal i = 0; i < size; ++i) {
		for (eslocal e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
			results[e * increment] = _value;
		}
	}
}


