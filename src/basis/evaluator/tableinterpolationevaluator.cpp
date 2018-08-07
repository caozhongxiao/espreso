
#include "tableinterpolationevaluator.h"

#include "../logging/logging.h"

using namespace espreso;

TableInterpolationEvaluator::TableInterpolationEvaluator(const std::vector<std::pair<double, double> > &table)
: _table(table)
{
	if (!table.size()) {
		ESINFO(GLOBAL_ERROR) << "Interpolation table with zero size.";
	}
}

void TableInterpolationEvaluator::evaluate(eslocal size, eslocal increment, const Point* cbegin, const double* tbegin, double time, double *results) const
{
	for (eslocal i = 0; i < size; ++i) {
		if (tbegin[i] < _table[0].first) {
			results[i * increment] = _table[0].second;
			break;
		}
		for (size_t j = 0; j < _table.size() - 1; j++) {
			if (_table[j].first < tbegin[i] && tbegin[i] < _table[j + 1].first) {
				double a = _table[j].first, b = _table[j + 1].first;
				double va = _table[j].second, vb = _table[j + 1].second;
				results[i * increment] = va + (vb - va) * (tbegin[i] - a) / (b - a);
				break;
			}
		}
		results[i * increment] = _table.back().second;
	}
}

std::string TableInterpolationEvaluator::getEXPRTKForm() const
{
	std::string exprtk;

	exprtk += "if((TEMPERATURE<" + std::to_string(_table.front().first) + "), " + std::to_string(_table.front().second) + ", 0) + ";

	for (size_t i = 0; i + 1 < _table.size(); i++) {
		exprtk += "if((TEMPERATURE<=" + std::to_string(_table[i].first) + " and " + std::to_string(_table[i + 1].first) + "<TEMPERATURE), ";
		double a = _table[i].first, b = _table[i + 1].first;
		double va = _table[i].second, vb = _table[i + 1].second;
		exprtk += std::to_string(va) + "+" + std::to_string(vb - va) + "*(TEMPERATURE - " + std::to_string(a) + "/" + std::to_string(b - a) + ")";
		exprtk += ", 0) + ";
	}

	exprtk += "if((" + std::to_string(_table.back().first) + "<TEMPERATURE), " + std::to_string(_table.back().second) + ", 0)";

	return exprtk;
}

