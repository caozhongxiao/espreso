
#ifndef SRC_BASIS_EVALUATORS_EVALUATOR_H_
#define SRC_BASIS_EVALUATORS_EVALUATOR_H_

#include "../expression/expression.h"

#include "omp.h"

#include <algorithm>
#include "../containers/point.h"

namespace espreso {

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
	static Evaluator* create(std::ifstream &is);

	virtual Evaluator* copy() const { return new Evaluator(*this); }

	virtual Type type() { return Type::DEFAULT; }
	virtual double evaluate(const Point &p, double time = 0, double temperature = 0, double pressure = 0, double velocity = 0) const { return 0; }

	virtual ~Evaluator() {};

	virtual bool isTimeDependent() const { return false; }
	virtual bool isTemperatureDependent() const { return false; }

	virtual void store(std::ofstream& os);
};

class ConstEvaluator: public Evaluator {

public:
	ConstEvaluator(double value);
	ConstEvaluator(std::ifstream &is);

	virtual Evaluator* copy() const { return new ConstEvaluator(*this); }

	virtual Type type() { return Type::CONST; }
	double inline evaluate(const Point &p, double time = 0, double temperature = 0, double pressure = 0, double velocity = 0) const { return _value; }

	virtual void store(std::ofstream& os);

private:
	double _value;
};


/**
 * Evaluator can be called from various threads.
 * Create one instance for each worker in order to avoid race conditions.
 */
class ExpressionEvaluator: public Evaluator {

public:
	ExpressionEvaluator(const std::string &expression, const std::vector<std::string> &variables = { "X", "Y", "Z", "TIME", "TEMPERATURE" });
	ExpressionEvaluator(std::ifstream &is);

	virtual Evaluator* copy() const { return new ExpressionEvaluator(*this); }

	virtual Type type() { return Type::EXPRESSION; }
	double inline evaluate(const Point &p, double time = 0, double temperature = 0, double pressure = 0, double velocity = 0) const
	{
		_values[omp_get_thread_num()][0] = p.x;
		_values[omp_get_thread_num()][1] = p.y;
		_values[omp_get_thread_num()][2] = p.z;
		_values[omp_get_thread_num()][3] = time;
		_values[omp_get_thread_num()][4] = temperature;
		_values[omp_get_thread_num()][5] = pressure;
		_values[omp_get_thread_num()][6] = velocity;
		return _expression[omp_get_thread_num()].evaluate(_values[omp_get_thread_num()]);
	}

	bool isTimeDependent() const { return _timeDependency; }
	bool isTemperatureDependent() const { return _temperatureDependency; }

	virtual void store(std::ofstream& os);

private:
	std::vector<Expression> _expression;
	std::vector<std::string> _variables;
	mutable std::vector<std::vector<double> >_values;
	bool _timeDependency, _temperatureDependency;
};

class TableEvaluator: public Evaluator {

public:
	enum class TableProperty {
		TIME,
		TEMPERATURE,
		PRESSURE,
		VELOCITY
	};

	TableEvaluator(
			const std::vector<std::vector<std::vector<double> > > &table,
			const std::vector<TableProperty> &properties,
			const std::vector<std::vector<double> > &axis);

	TableEvaluator(std::ifstream &is);

	virtual Evaluator* copy() const { return new TableEvaluator(*this); }

	virtual Type type() { return Type::TABLE; }
	virtual double inline evaluate(const Point &p, double time = 0, double temperature = 0, double pressure = 0, double velocity = 0) const
	{
		std::vector<size_t> cell(_dimension);

		for (size_t i = 0; i < _dimension; i++) {
			switch (_properties[i]) {
			case TableProperty::TIME:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), time) - _axis[i].begin();
				break;
			case TableProperty::TEMPERATURE:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), temperature) - _axis[i].begin();
				break;
			case TableProperty::PRESSURE:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), pressure) - _axis[i].begin();
				break;
			case TableProperty::VELOCITY:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), velocity) - _axis[i].begin();
				break;
			}
		}
		return _table[_dimension > 0 ? cell[0] : 0][_dimension > 1 ? cell[1] : 0][_dimension > 2 ? cell[2] : 0];
	}

	bool isTimeDependent() const { return std::any_of(_properties.begin(), _properties.end(), [] (const TableProperty &p) { return p == TableProperty::TIME; }); }
	bool isTemperatureDependent() const { return std::any_of(_properties.begin(), _properties.end(), [] (const TableProperty &p) { return p == TableProperty::TEMPERATURE; }); }

	virtual void store(std::ofstream& os);

protected:
	size_t _dimension;
	std::vector<std::vector<std::vector<double> > > _table;
	std::vector<TableProperty> _properties;
	std::vector<std::vector<double> > _axis;
};

class TableInterpolationEvaluator: public Evaluator {

public:
	TableInterpolationEvaluator(const std::vector<std::pair<double, double> > &table);
	TableInterpolationEvaluator(std::ifstream &is);

	virtual Evaluator* copy() const { return new TableInterpolationEvaluator(*this); }

	virtual Type type() { return Type::TABLE_INTERPOLATION; }
	virtual double inline evaluate(const Point &p, double time = 0, double temperature = 0, double pressure = 0, double velocity = 0) const
	{
		if (temperature < table[0].first) {
			return table[0].second;
		}
		for (size_t i = 0; i < table.size() - 1; i++) {
			if (table[i].first < temperature && temperature < table[i + 1].first) {
				double a = table[i].first  , b = table[i + 1].first;
				double va = table[i].second, vb = table[i + 1].second;
				return va + (vb - va) * (temperature - a) / (b - a);
			}
		}
		return table.back().second;
	}

	bool isTimeDependent() const { return false; }
	bool isTemperatureDependent() const { return true; }

	virtual void store(std::ofstream& os);

	std::vector<std::pair<double, double> > table;
};

namespace input { class API; }

class ArrayEvaluator: public Evaluator {

	friend class input::API;
public:
	ArrayEvaluator(std::vector<eslocal> &indices, std::vector<double> &values, eslocal offset);
	ArrayEvaluator(eslocal size, eslocal *indices, double *values, eslocal offset);

	void addIndex(eslocal index, eslocal value);

	virtual Evaluator* copy() const { return new ArrayEvaluator(*this); }

	virtual Type type() { return Type::ARRAY; }
	virtual double evaluate(const Point &p, double time = 0, double temperature = 0, double pressure = 0, double velocity = 0) const;
	virtual double evaluate(eslocal index) const;

	virtual void store(std::ofstream& os);

protected:
	std::vector<eslocal> _indices;
	std::vector<double> _values;
};

}



#endif /* SRC_BASIS_EVALUATORS_EVALUATOR_H_ */
