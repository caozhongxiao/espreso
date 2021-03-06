
#include "tableinterpolationevaluator.h"

#include "esinfo/eslog.hpp"

using namespace espreso;

TableInterpolationEvaluator::TableInterpolationEvaluator(const std::vector<std::pair<double, double> > &table)
: _table(table)
{
	if (!table.size()) {
		eslog::globalerror("Interpolation table with zero size.\n");
	}
}

void TableInterpolationEvaluator::evalVector(esint size, esint increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	for (esint i = 0; i < size; ++i) {
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

void TableInterpolationEvaluator::evalFiltered(esint size, esint increment, esint *elements, esint *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		for (esint e = distribution[elements[i]]; e < distribution[elements[i] + 1]; ++e) {
			if (tbegin[e] < _table[0].first) {
				results[e] = _table[0].second;
				break;
			}
			for (size_t j = 0; j < _table.size() - 1; j++) {
				if (_table[j].first < tbegin[i] && tbegin[i] < _table[j + 1].first) {
					double a = _table[j].first, b = _table[j + 1].first;
					double va = _table[j].second, vb = _table[j + 1].second;
					results[e] = va + (vb - va) * (tbegin[i] - a) / (b - a);
					break;
				}
			}
			results[e * increment] = _table.back().second;
		}
	}
}

void TableInterpolationEvaluator::evalSelected(esint size, esint increment, esint *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const
{
	for (esint i = 0; i < size; ++i) {
		if (tbegin[selection[i]] < _table[0].first) {
			results[i * increment] = _table[0].second;
			break;
		}
		for (size_t j = 0; j < _table.size() - 1; j++) {
			if (_table[j].first < tbegin[selection[i]] && tbegin[selection[i]] < _table[j + 1].first) {
				double a = _table[j].first, b = _table[j + 1].first;
				double va = _table[j].second, vb = _table[j + 1].second;
				results[i * increment] = va + (vb - va) * (tbegin[selection[i]] - a) / (b - a);
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

