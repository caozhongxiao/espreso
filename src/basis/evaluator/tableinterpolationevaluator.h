
#ifndef SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_

#include "evaluator.h"

#include <utility>
#include <vector>

namespace espreso {

class TableInterpolationEvaluator: public Evaluator {

public:
	TableInterpolationEvaluator(const std::vector<std::pair<double, double> > &table);

	Type type() { return Type::TABLE_INTERPOLATION; }
	Evaluator* copy() const { return new TableInterpolationEvaluator(*this); }

	void evaluate(eslocal size, eslocal increment, const Point* cbegin, const double* tbegin, double time, double *results) const;

	bool isCoordinateDependent() const { return false; }
	bool isTimeDependent() const { return false; }
	bool isTemperatureDependent() const { return true; }

protected:
	std::vector<std::pair<double, double> > _table;
};

}



#endif /* SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_ */
