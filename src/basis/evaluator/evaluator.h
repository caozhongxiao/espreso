
#ifndef SRC_BASIS_EVALUATOR_EVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EVALUATOR_H_

#include <string>

namespace espreso {

struct Point;

enum EvaluatorParameters: int {
	VALUE       = 0, // never const
	COORDINATE  = 1 << 0,
	TIME        = 1 << 1,
	TEMPERATURE = 1 << 2
};

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

	virtual void evalVector(esint size, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		evalVector(size, 1, csize, cbegin, tbegin, time, results);
	}
	virtual void evalVector(esint size, esint increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		for (esint i = 0; i < size; ++i) {
			results[i * increment] = 0;
		}
	};

	virtual void evalFiltered(esint size, esint *elements, esint *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		evalFiltered(size, 1, elements, distribution, csize, cbegin, tbegin, time, results);
	}

	virtual void evalFiltered(esint size, esint increment, esint *elements, esint *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		for (esint i = 0; i < size; ++i) {
			for (esint e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
				results[e * increment] = 0;
			}
		}
	}

	virtual void evalSelected(esint size, esint *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		evalSelected(size, 1, selection, csize, cbegin, tbegin, time, results);
	}

	virtual void evalSelected(esint size, esint increment, esint *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
	{
		for (esint i = 0; i < size; ++i) {
			results[i * increment] = 0;
		}
	}

	virtual double evaluate(double r) const { return 0; }

	// TODO: remove
	double evaluate(const Point &p, double time, double temperature) const;

	virtual bool isConstant() const { return true; }
	virtual bool isCoordinateDependent() const { return false; }
	virtual bool isTimeDependent() const { return false; }
	virtual bool isTemperatureDependent() const { return false; }

	virtual bool isConstant(EvaluatorParameters parameters) const
	{
		return
				(parameters != EvaluatorParameters::VALUE) &&
				(!(parameters & EvaluatorParameters::COORDINATE) || !isCoordinateDependent()) &&
				(!(parameters & EvaluatorParameters::TIME) || !isTimeDependent()) &&
				(!(parameters & EvaluatorParameters::TEMPERATURE) || !isTemperatureDependent());
	}

	virtual std::string getEXPRTKForm() const { return "0"; }
};

}



#endif /* SRC_BASIS_EVALUATOR_EVALUATOR_H_ */
