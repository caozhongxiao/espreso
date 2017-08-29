
#include <fstream>

#include "evaluator.h"
#include <numeric>
#include "../../config/ecf/environment.h"
#include "../../basis/utilities/parser.h"

using namespace espreso;


Evaluator* Evaluator::create(std::ifstream &is)
{
	Evaluator::Type type;
	is.read(reinterpret_cast<char *>(&type), sizeof(Evaluator::Type));

	switch (type) {
	case Type::DEFAULT: return new Evaluator();
	case Type::CONST: return new ConstEvaluator(is);
	case Type::EXPRESSION: return new ExpressionEvaluator(is);
	case Type::TABLE: return new TableEvaluator(is);
	case Type::TABLE_INTERPOLATION: return new TableInterpolationEvaluator(is);
	case Type::ARRAY: ESINFO(GLOBAL_ERROR) << "Implement loading of Array evaluator"; return NULL;
	default: ESINFO(GLOBAL_ERROR) << "Unknown evaluator type"; return NULL;
	}
}

void Evaluator::store(std::ofstream& os)
{
	Type type = Type::DEFAULT;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
}

ConstEvaluator::ConstEvaluator(double value)
: _value(value)
{

}

ConstEvaluator::ConstEvaluator(std::ifstream &is)
{
	is.read(reinterpret_cast<char *>(&_value), sizeof(double));
}

void ConstEvaluator::store(std::ofstream& os)
{
	Type type = Type::CONST;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
	os.write(reinterpret_cast<const char *>(&_value), sizeof(double));
}

ExpressionEvaluator::ExpressionEvaluator(const std::string &expression, const std::vector<std::string> &variables)
: _variables({ "X", "Y", "Z", "TIME", "TEMPERATURE", "PRESSURE", "VELOCITY" })
{
	for (size_t i = 0; i < variables.size(); i++) {
		if (std::find(_variables.begin(), _variables.end(), variables[i]) == _variables.end()) {
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: ExpressionEvaluator not supports variable: '" << variables[i] << "'.";
		}
	}
	_expression.resize(environment->OMP_NUM_THREADS, Expression(expression, _variables));
	_values.resize(environment->OMP_NUM_THREADS, std::vector<double>(_variables.size()));
	_timeDependency = StringCompare::contains(expression, { "TIME" });
	_temperatureDependency = StringCompare::contains(expression, { "TEMPERATURE" });
}

ExpressionEvaluator::ExpressionEvaluator(std::ifstream &is)
{
	eslocal size;
	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	char *buffer = new char[size];
	is.read(buffer, size);
	std::string expression(buffer, size);
	delete buffer;
	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
	for (eslocal i = 0; i < size; i++) {
		eslocal vsize;
		is.read(reinterpret_cast<char *>(&vsize), sizeof(eslocal));
		buffer = new char[vsize];
		is.read(buffer, size);
		_variables.push_back(std::string(buffer, vsize));
		delete buffer;
	}
	_expression.resize(environment->OMP_NUM_THREADS, Expression(expression, _variables));
	_values.resize(environment->OMP_NUM_THREADS, std::vector<double>(_variables.size()));
	_timeDependency = StringCompare::contains(_expression[0].expression(), { "TIME" });
	_temperatureDependency = StringCompare::contains(_expression[0].expression(), { "TEMPERATURE" });
}

void ExpressionEvaluator::store(std::ofstream& os)
{
	Type type = Type::EXPRESSION;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
	eslocal size = _expression[0].expression().size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(eslocal));
	os.write(_expression[0].expression().c_str(), _expression[0].expression().size());
	size = _variables.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(eslocal));
	for (size_t i = 0; i < _variables.size(); i++) {
		size = _variables[i].size();
		os.write(reinterpret_cast<const char *>(&size), sizeof(eslocal));
		os.write(_variables[i].c_str(), _variables[i].size());
	}
}


TableEvaluator::TableEvaluator(
			const std::vector<std::vector<std::vector<double> > > &table,
			const std::vector<TableProperty> &properties,
			const std::vector<std::vector<double> > &axis)
: _dimension(properties.size()), _table(table), _properties(properties), _axis(axis)
{

}

TableEvaluator::TableEvaluator(std::ifstream &is)
{
	is.read(reinterpret_cast<char *>(&_dimension), sizeof(size_t));

	size_t size;
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	_table.resize(size);
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _table.size(); i++) {
		_table[i].resize(size);
	}
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _table.size(); i++) {
		for (size_t j = 0; j < _table[i].size(); j++) {
			_table[i][j].resize(size);
			for (size_t k = 0; k < _table[i][j].size(); k++) {
				is.read(reinterpret_cast<char *>(&_table[i][j][k]), sizeof(double));
			}
		}
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	_properties.resize(size);
	for (size_t i = 0; i < _properties.size(); i++) {
		is.read(reinterpret_cast<char *>(&_properties[i]), sizeof(TableProperty));
	}

	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	_axis.resize(size);
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _axis.size(); i++) {
		_axis[i].resize(size);
		for (size_t j = 0; j < _axis[i].size(); j++) {
			is.read(reinterpret_cast<char *>(&_axis[i][j]), sizeof(double));
		}
	}
}

void TableEvaluator::store(std::ofstream& os)
{
	Type type = Type::TABLE;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));
	size_t size = _dimension;
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));

	size = _table.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	size = _table[0].size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	size = _table[0][0].size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _table.size(); i++) {
		for (size_t j = 0; j < _table[i].size(); j++) {
			for (size_t k = 0; k < _table[i][j].size(); k++) {
				os.write(reinterpret_cast<const char *>(&_table[i][j][k]), sizeof(double));
			}
		}
	}

	size = _properties.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _properties.size(); i++) {
		os.write(reinterpret_cast<const char *>(&_properties[i]), sizeof(TableProperty));
	}

	size = _axis.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	size = _axis[0].size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	for (size_t i = 0; i < _axis.size(); i++) {
		for (size_t j = 0; j < _axis[i].size(); j++) {
			os.write(reinterpret_cast<const char *>(&_axis[i][j]), sizeof(double));
		}
	}
}

TableInterpolationEvaluator::TableInterpolationEvaluator(const std::vector<std::pair<double, double> > &table)
: table(table)
{
	if (!table.size()) {
		ESINFO(GLOBAL_ERROR) << "Interpolation table with zero size.";
	}
}

TableInterpolationEvaluator::TableInterpolationEvaluator(std::ifstream &is)
{
	size_t size;
	is.read(reinterpret_cast<char *>(&size), sizeof(size_t));
	table.resize(size);
	is.read(reinterpret_cast<char *>(table.data()), 2 * sizeof(size_t) * table.size());
}

void TableInterpolationEvaluator::store(std::ofstream& os)
{
	Type type = Type::TABLE_INTERPOLATION;
	os.write(reinterpret_cast<const char *>(&type), sizeof(Evaluator::Type));

	size_t size = table.size();
	os.write(reinterpret_cast<const char *>(&size), sizeof(size_t));
	os.write(reinterpret_cast<const char *>(table.data()), 2 * sizeof(double) * table.size());
}

ArrayEvaluator::ArrayEvaluator(std::vector<eslocal> &indices, std::vector<double> &values, eslocal offset)
{
	std::vector<eslocal> permutation(indices.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return indices[i] < indices[j]; });
	_indices.reserve(indices.size());
	_values.reserve(indices.size());
	for (size_t i = 0; i < permutation.size(); i++) {
		_indices.push_back(indices[permutation[i]] - offset);
		_values.push_back(values[permutation[i]]);
	}
}

ArrayEvaluator::ArrayEvaluator(eslocal size, eslocal *indices, double *values, eslocal offset)
{
	std::vector<eslocal> permutation(size);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return indices[i] < indices[j]; });
	_indices.reserve(size);
	_values.reserve(size);
	for (size_t i = 0; i < permutation.size(); i++) {
		_indices.push_back(indices[permutation[i]] - offset);
		_values.push_back(values[permutation[i]]);
	}
}

void ArrayEvaluator::addIndex(eslocal index, eslocal value)
{
	size_t offset = std::lower_bound(_indices.begin(), _indices.end(), index) - _indices.begin();
	_indices.insert(_indices.begin() + offset, index);
	_values.insert(_values.begin() + offset, value);
}


void ArrayEvaluator::store(std::ofstream& os)
{
	ESINFO(GLOBAL_ERROR) << "Implement store ArrayEvaluator.";
}

double inline ArrayEvaluator::evaluate(const Point &p, double time, double temperature, double pressure, double velocity) const
{
	ESINFO(ERROR) << "Invalid calling of ArrayEvaluator.";
	return 0;
}

double inline ArrayEvaluator::evaluate(eslocal index) const
{
	auto it = std::lower_bound(_indices.begin(), _indices.end(), index);
	if (it != _indices.end() && *it == index) {
		return _values[it - _indices.begin()];
	} else {
		ESINFO(ERROR) << "Array evaluator has no specified value for index '" << index << "'";
		return 0;
	}
}




