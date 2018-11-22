
#include "tableevaluator.h"

#include <algorithm>
#include <functional>

using namespace espreso;


TableEvaluator::TableEvaluator(
		const std::vector<std::vector<std::vector<double> > > &table,
		const std::vector<TableProperty> &properties,
		const std::vector<std::vector<double> > &axis)
: _dimension(properties.size()), _table(table), _properties(properties), _axis(axis)
{
	_temperatureDependency = std::any_of(_properties.begin(), _properties.end(), [] (const TableProperty &p) { return p == TableProperty::TEMPERATURE; });
	_timeDependency = std::any_of(_properties.begin(), _properties.end(), [] (const TableProperty &p) { return p == TableProperty::TIME; });
}

void TableEvaluator::evalVector(eslocal size, eslocal increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	std::vector<size_t> cell(_dimension);

	for (eslocal j = 0; j < size; ++j) {
		for (size_t i = 0; i < _dimension; i++) {
			switch (_properties[i]) {
			case TableProperty::TEMPERATURE:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), tbegin[j]) - _axis[i].begin();
				break;
			case TableProperty::TIME:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), time) - _axis[i].begin();
				break;

			}
		}
		results[j * increment] = _table[_dimension > 0 ? cell[0] : 0][_dimension > 1 ? cell[1] : 0][_dimension > 2 ? cell[2] : 0];
	}
}

void TableEvaluator::evalFiltered(eslocal size, eslocal increment, eslocal *elements, eslocal *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	std::vector<size_t> cell(_dimension);

	for (eslocal j = 0; j < size; ++j) {
		for (eslocal e = distribution[elements[j]]; e < distribution[elements[j] + 1]; ++e) {
			for (size_t i = 0; i < _dimension; i++) {
				switch (_properties[i]) {
				case TableProperty::TEMPERATURE:
					cell[i] = std::find(_axis[i].begin(), _axis[i].end(), tbegin[e]) - _axis[i].begin();
					break;
				case TableProperty::TIME:
					cell[i] = std::find(_axis[i].begin(), _axis[i].end(), time) - _axis[i].begin();
					break;

				}
			}
			results[e * increment] = _table[_dimension > 0 ? cell[0] : 0][_dimension > 1 ? cell[1] : 0][_dimension > 2 ? cell[2] : 0];
		}
	}
}

void TableEvaluator::evalSelected(eslocal size, eslocal increment, eslocal *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	std::vector<size_t> cell(_dimension);

	for (eslocal j = 0; j < size; ++j) {
		for (size_t i = 0; i < _dimension; i++) {
			switch (_properties[i]) {
			case TableProperty::TEMPERATURE:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), tbegin[selection[j]]) - _axis[i].begin();
				break;
			case TableProperty::TIME:
				cell[i] = std::find(_axis[i].begin(), _axis[i].end(), time) - _axis[i].begin();
				break;

			}
		}
		results[j * increment] = _table[_dimension > 0 ? cell[0] : 0][_dimension > 1 ? cell[1] : 0][_dimension > 2 ? cell[2] : 0];
	}
}

std::string TableEvaluator::getEXPRTKForm() const
{
	std::string exprtk;

	auto getproperty = [&] (TableProperty tproperty) {
		switch (tproperty) {
		case TableProperty::TEMPERATURE:
			return "TEMPERATURE";
		case TableProperty::TIME:
			return "TIME";
		default:
			return "";
		}
	};

	std::vector<size_t> cell(_dimension);

	auto eval = [&] () {
		return _table[_dimension > 0 ? cell[0] : 0][_dimension > 1 ? cell[1] : 0][_dimension > 2 ? cell[2] : 0];
	};

	std::function<void(size_t, const std::string&, const std::vector<double>&)> addAxis = [&] (size_t dimension, const std::string &property, const std::vector<double> &values) {
		exprtk += "if((" + property + "<" + std::to_string(values.front()) + "), ";
		if (dimension + 1 == _dimension) {
			exprtk += std::to_string(eval()) + ", 0";
		} else {
			addAxis(dimension + 1, getproperty(_properties[dimension + 1]), _axis[dimension + 1]);
		}
		exprtk += ") + ";

		for (size_t i = 0; i + 1 < values.size(); i++) {
			exprtk += "if((" + property + "<=" + std::to_string(values[i]) + " and " + property + "<" + std::to_string(values[i + 1]) + "), ";
			if (dimension + 1 == _dimension) {
				exprtk += std::to_string(eval()) + ", 0";
			} else {
				addAxis(dimension + 1, getproperty(_properties[dimension + 1]), _axis[dimension + 1]);
			}
			exprtk += ") + ";
		}

		exprtk += "if((" + std::to_string(values.back()) + "<" + property + "), ";
		if (dimension + 1 == _dimension) {
			exprtk += std::to_string(eval()) + ", 0";
		} else {
			addAxis(dimension + 1, getproperty(_properties[dimension + 1]), _axis[dimension + 1]);
		}
		exprtk += ")";
	};

	addAxis(0, getproperty(_properties[0]), _axis[0]);

	return exprtk;
}

