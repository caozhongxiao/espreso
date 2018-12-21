
#ifndef SRC_BASIS_EVALUATOR_TABLEEVALUATOR_H_
#define SRC_BASIS_EVALUATOR_TABLEEVALUATOR_H_

#include "evaluator.h"

#include <cstddef>
#include <vector>

namespace espreso {

class TableEvaluator: public Evaluator {

public:
	enum class TableProperty {
		TEMPERATURE,
		TIME
	};

	TableEvaluator(
			const std::vector<std::vector<std::vector<double> > > &table,
			const std::vector<TableProperty> &properties,
			const std::vector<std::vector<double> > &axis);

	virtual Type type() { return Type::TABLE; }
	virtual Evaluator* copy() const { return new TableEvaluator(*this); }

	void evalVector(esint size, esint increment, int csize, const double* cbegin, const double* tbegin, double time, double *results) const;
	void evalFiltered(esint size, esint increment, esint *elements, esint *distribution, int csize, const double* cbegin, const double* tbegin, double time, double *results) const;
	void evalSelected(esint size, esint increment, esint *selection, int csize, const double* cbegin, const double* tbegin, double time, double *results) const;

	bool isConstant() const { return !_timeDependency && !_temperatureDependency; }
	bool isCoordinateDependent() const { return false; }
	bool isTimeDependent() const { return _timeDependency; }
	bool isTemperatureDependent() const { return _temperatureDependency; }

	std::string getEXPRTKForm() const;

protected:
	size_t _dimension;
	std::vector<std::vector<std::vector<double> > > _table;
	std::vector<TableProperty> _properties;
	std::vector<std::vector<double> > _axis;
	bool _temperatureDependency, _timeDependency;
};

}

#endif /* SRC_BASIS_EVALUATOR_TABLEEVALUATOR_H_ */
