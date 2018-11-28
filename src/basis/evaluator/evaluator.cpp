
#include "evaluator.h"

#include "../containers/point.h"

using namespace espreso;

double Evaluator::evaluate(const Point &p, double time, double temperature) const
{
	double result;
	evalVector(1, 3, &p.x, &temperature, time, &result);
	return result;
}



