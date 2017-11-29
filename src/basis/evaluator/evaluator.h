
#ifndef SRC_BASIS_EVALUATOR_EVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EVALUATOR_H_

namespace espreso {

struct Point;

class Evaluator {

protected:
	enum class Type: int {
		DEFAULT,
		CONST,
		EXPRESSION,
		TABLE,
		TABLE_INTERPOLATION,
		ARRAY
	};

public:
	virtual Type type() { return Type::DEFAULT; }
	virtual Evaluator* copy() const { return new Evaluator(*this); }

	virtual ~Evaluator() {};

	virtual void evaluate(eslocal size, const Point* cbegin, const double* tbegin, double time, double *results) const {};

	// TODO: remove
	double evaluate(const Point &p, double temperature, double time);

	virtual bool isCoordinateDependent() const { return false; }
	virtual bool isTimeDependent() const { return false; }
	virtual bool isTemperatureDependent() const { return false; }
};

}



#endif /* SRC_BASIS_EVALUATOR_EVALUATOR_H_ */
