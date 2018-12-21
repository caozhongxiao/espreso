
#include "constevaluator.h"

using namespace espreso;

void ConstEvaluator::evalVector(esint size, esint increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		results[i * increment] = _value;
	}
}


void ConstEvaluator::evalFiltered(esint size, esint increment, esint *elements, esint *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		for (esint e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
			results[e * increment] = _value;
		}
	}
}


