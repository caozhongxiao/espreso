
#ifndef SRC_BASIS_EVALUATOR_EVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EVALUATOR_H_

#include <string>

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

	virtual void evalVector(eslocal size, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		evalVector(size, 1, csize, cbegin, tbegin, time, results);
	}
	virtual void evalVector(eslocal size, eslocal increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		for (eslocal i = 0; i < size; ++i) {
			results[i * increment] = 0;
		}
	};

	virtual void evalFiltered(eslocal size, eslocal *elements, eslocal *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		evalFiltered(size, 1, elements, distribution, csize, cbegin, tbegin, time, results);
	}

	virtual void evalFiltered(eslocal size, eslocal increment, eslocal *elements, eslocal *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		for (eslocal i = 0; i < size; ++i) {
			for (eslocal e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
				results[e * increment] = 0;
			}
		}
	}

	virtual void evalSelected(eslocal size, eslocal *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		evalSelected(size, 1, selection, csize, cbegin, tbegin, time, results);
	}

	virtual void evalSelected(eslocal size, eslocal increment, eslocal *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		for (eslocal i = 0; i < size; ++i) {
			results[i] = 0;
		}
	}

	virtual double evaluate(double r) const { return 0; }

	// TODO: remove
	double evaluate(const Point &p, double time, double temperature) const;

	virtual bool isConstant() const { return true; }
	virtual bool isCoordinateDependent() const { return false; }
	virtual bool isTimeDependent() const { return false; }
	virtual bool isTemperatureDependent() const { return false; }

	virtual std::string getEXPRTKForm() const { return "0"; }
};

}



#endif /* SRC_BASIS_EVALUATOR_EVALUATOR_H_ */
