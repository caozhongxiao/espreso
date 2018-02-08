
#include "expressionevaluator.h"
#include "../expression/expression.h"

#include "omp.h"

#include "../../basis/containers/point.h"
#include "../../basis/utilities/parser.h"
#include "../../config/ecf/environment.h"

using namespace espreso;

ExpressionEvaluator::ExpressionEvaluator(const std::string &expression)
{
	size_t threads = environment->OMP_NUM_THREADS;

	_expressions.resize(environment->OMP_NUM_THREADS);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		_expressions[t] = new Expression(expression, ExpressionEvaluator::variables());
	}

	_coordinateDependency = StringCompare::contains(expression, { "X", "Y", "Z" });
	_timeDependency = StringCompare::contains(expression, { "TIME" });
	_temperatureDependency = StringCompare::contains(expression, { "TEMPERATURE" });
}

ExpressionEvaluator::ExpressionEvaluator(const ExpressionEvaluator &other)
{
	size_t threads = environment->OMP_NUM_THREADS;

	_expressions.resize(environment->OMP_NUM_THREADS);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		_expressions[t] = new Expression(*other._expressions[t]);
	}

	_coordinateDependency = other._coordinateDependency;
	_timeDependency = other._timeDependency;
	_temperatureDependency = other._temperatureDependency;
}

void ExpressionEvaluator::evaluate(eslocal size, eslocal increment, const Point* cbegin, const double* tbegin, double time, double *results) const
{
	int thread = omp_get_thread_num();
	for (eslocal i = 0; i < size; ++i) {
		if (cbegin != NULL) {
			_expressions[thread]->values[0] = cbegin[i].x;
			_expressions[thread]->values[1] = cbegin[i].y;
			_expressions[thread]->values[2] = cbegin[i].z;
		}
		if (tbegin != NULL) {
			_expressions[thread]->values[3] = tbegin[i];
		}
		_expressions[thread]->values[4] = time;
		results[i * increment] = _expressions[thread]->evaluate();
	}
}

double ExpressionEvaluator::evaluate(double r) const
{
	int thread = omp_get_thread_num();
	_expressions[thread]->values[5] = r;
	return _expressions[thread]->evaluate();
}

std::string ExpressionEvaluator::getEXPRTKForm() const
{
	return _expressions.front()->expression();
}




