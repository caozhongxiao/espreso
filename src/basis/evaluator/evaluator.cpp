
#include "evaluator.h"

using namespace espreso;


double Evaluator::evaluate(const Point &p, double temperature, double time)
{
	double result;
	evaluate(1, &p, &temperature, time, &result);
	return result;
}



