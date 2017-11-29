
#include "tableevaluator.h"

#include <algorithm>

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

void TableEvaluator::evaluate(eslocal size, const Point* cbegin, const double* tbegin, double time, double *results) const
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
		results[j] = _table[_dimension > 0 ? cell[0] : 0][_dimension > 1 ? cell[1] : 0][_dimension > 2 ? cell[2] : 0];
	}
}

