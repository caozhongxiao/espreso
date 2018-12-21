
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

	void evalVector(esint size, esint increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const;
	void evalFiltered(esint size, esint increment, esint *elements, esint *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const;
	void evalSelected(esint size, esint increment, esint *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const;

	bool isConstant() const { return false; }
	bool isCoordinateDependent() const { return false; }
	bool isTimeDependent() const { return false; }
	bool isTemperatureDependent() const { return true; }

	std::string getEXPRTKForm() const;

protected:
	std::vector<std::pair<double, double> > _table;
};

}



#endif /* SRC_BASIS_EVALUATOR_TABLEINTERPOLATIONEVALUATOR_H_ */
