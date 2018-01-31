
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

	virtual void evaluate(eslocal size, const Point* cbegin, const double* tbegin, double time, double *results) const
	{
		evaluate(size, 1, cbegin, tbegin, time, results);
	}
	virtual void evaluate(eslocal size, eslocal increment, const Point* cbegin, const double* tbegin, double time, double *results) const {};

	virtual double evaluate(double r) const { return 0; }

	// TODO: remove
	double evaluate(const Point &p, double time, double temperature) const;

	virtual bool isCoordinateDependent() const { return false; }
	virtual bool isTimeDependent() const { return false; }
	virtual bool isTemperatureDependent() const { return false; }
};

}



#endif /* SRC_BASIS_EVALUATOR_EVALUATOR_H_ */
