
#include "esinfo/envinfo.h"
#include "expressionevaluator.h"

#include "basis/containers/point.h"
#include "basis/utilities/parser.h"
#include "basis/expression/expression.h"
#include "omp.h"

using namespace espreso;

ExpressionEvaluator::ExpressionEvaluator(const std::string &expression)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	_expressions.resize(info::env::OMP_NUM_THREADS);

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
	size_t threads = info::env::OMP_NUM_THREADS;

	_expressions.resize(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		_expressions[t] = new Expression(*other._expressions[t]);
	}

	_coordinateDependency = other._coordinateDependency;
	_timeDependency = other._timeDependency;
	_temperatureDependency = other._temperatureDependency;
}

void ExpressionEvaluator::evalVector(esint size, esint increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		if (cbegin != NULL) {
			_expressions[thread]->values[0] = cbegin[i * csize + 0];
			_expressions[thread]->values[1] = cbegin[i * csize + 1];
			_expressions[thread]->values[2] = csize == 3 ? cbegin[i * csize + 2] : 0;
		}
		if (tbegin != NULL) {
			_expressions[thread]->values[3] = tbegin[i];
		}
		_expressions[thread]->values[4] = time;
		results[i * increment] = _expressions[thread]->evaluate();
	}
}

void ExpressionEvaluator::evalFiltered(esint size, esint increment, esint *elements, esint *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		for (esint e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
			if (cbegin != NULL) {
				_expressions[thread]->values[0] = cbegin[e * csize + 0];
				_expressions[thread]->values[1] = cbegin[e * csize + 1];
				_expressions[thread]->values[2] = csize == 3 ? cbegin[e * csize + 2] : 0;
			}
			if (tbegin != NULL) {
				_expressions[thread]->values[3] = tbegin[e];
			}
			_expressions[thread]->values[4] = time;
			results[e * increment] = _expressions[thread]->evaluate();
		}
	}
}

void ExpressionEvaluator::evalSelected(esint size, esint increment, esint *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	int thread = omp_get_thread_num();
	for (esint i = 0; i < size; ++i) {
		if (cbegin != NULL) {
			_expressions[thread]->values[0] = cbegin[selection[i] * csize + 0];
			_expressions[thread]->values[1] = cbegin[selection[i] * csize + 1];
			_expressions[thread]->values[2] = csize == 3 ? cbegin[selection[i] * csize + 2] : 0;
		}
		if (tbegin != NULL) {
			_expressions[thread]->values[3] = tbegin[selection[i]];
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




