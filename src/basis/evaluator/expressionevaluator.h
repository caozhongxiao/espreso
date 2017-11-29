
#ifndef SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_

#include "evaluator.h"

#include <string>
#include <vector>

namespace espreso {

class Expression;

/**
 * Evaluator can be called from various threads.
 * Create one instance for each worker in order to avoid race conditions.
 */
class ExpressionEvaluator: public Evaluator {

public:
	ExpressionEvaluator(const std::string &expression);
	ExpressionEvaluator(const ExpressionEvaluator &other);

	Type type() { return Type::EXPRESSION; }
	virtual Evaluator* copy() const { return new ExpressionEvaluator(*this); }

	void evaluate(eslocal size, const Point* cbegin, const double* tbegin, double time, double *results) const;

	bool isCoordinateDependent() const { return _coordinateDependency; }
	bool isTimeDependent() const { return _timeDependency; }
	bool isTemperatureDependent() const { return _temperatureDependency; }

protected:


	std::vector<Expression*> _expressions;
	bool _coordinateDependency, _temperatureDependency, _timeDependency;
};

}


#endif /* SRC_BASIS_EVALUATOR_EXPRESSIONEVALUATOR_H_ */
