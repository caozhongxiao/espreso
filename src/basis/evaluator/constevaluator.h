
#ifndef SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_

#include "evaluator.h"

namespace espreso {

class ConstEvaluator: public Evaluator {

public:
	ConstEvaluator(double value): _value(value) {}

	Type type() { return Type::CONST; }
	virtual Evaluator* copy() const { return new ConstEvaluator(*this); }

	void evaluate(eslocal size, const Point* cbegin, const double* tbegin, double time, double *results) const;

protected:
	double _value;
};

}


#endif /* SRC_BASIS_EVALUATOR_CONSTEVALUATOR_H_ */
