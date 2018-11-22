
#ifndef SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_

#include "evaluator.h"

namespace espreso {

class ConstEvaluator: public Evaluator {

public:
	ConstEvaluator(double value): _value(value) {}

	Type type() { return Type::CONST; }
	virtual Evaluator* copy() const { return new ConstEvaluator(*this); }

	void evalVector(eslocal size, eslocal increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const;
	void evalFiltered(eslocal size, eslocal increment, eslocal *elements, eslocal *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const;

	void evalSelected(eslocal size, eslocal increment, eslocal *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		evalVector(size, increment, csize, cbegin, tbegin, time, results);
	}

	double evaluate(double r) const { return _value; }

	std::string getEXPRTKForm() const { return std::to_string(_value); }

protected:
	double _value;
};

}


#endif /* SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_ */
