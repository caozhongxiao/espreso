
#include "controler.h"

#include "../../../config/ecf/environment.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../config/expression.h"
#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"

using namespace espreso;


Controler::Controler(Mesh &mesh, const Step &step)
: _mesh(mesh), _step(step)
{
	size_t threads = environment->OMP_NUM_THREADS;
	_nDistribution.resize(threads + 1);

	for (size_t t = 0; t < threads; t++) {
		_nDistribution[t + 1] = _nDistribution[t];
		for (eslocal d = _mesh.elements->domainDistribution[t]; d != _mesh.elements->domainDistribution[t + 1]; ++d) {
			eslocal dbegin = _mesh.elements->elementsDistribution[d];
			eslocal dend = _mesh.elements->elementsDistribution[d + 1];
			_nDistribution[t + 1] += (_mesh.elements->procNodes->begin() + dend)->begin() - (_mesh.elements->procNodes->begin() + dbegin)->begin();
		}
	}
}

Controler::Parameter::~Parameter()
{
	if (data != NULL) {
		delete data;
	}
}

bool Controler::tryElementConstness(const std::map<std::string, ECFExpression> &values, double &defaultValue)
{
	if (values.size() == 0) {
		return true;
	}
	if (values.size() > 1) {
		return false;
	}
	if (_mesh.onAllElements(values.begin()->first) && values.begin()->second.evaluator->isConstant()) {
		values.begin()->second.evaluator->evaluate(1, 0, NULL, NULL, 0, &defaultValue);
		return true;
	}
	return false;
}

bool Controler::tryElementConstness(const std::map<std::string, ECFExpressionVector> &values, Point &defaultValue)
{
	if (values.size() == 0) {
		return true;
	}
	if (values.size() > 1) {
		return false;
	}
	if (_mesh.onAllElements(values.begin()->first)) {
		if (
				values.begin()->second.x.evaluator->isConstant() &&
				values.begin()->second.y.evaluator->isConstant() &&
				values.begin()->second.z.evaluator->isConstant()) {

			values.begin()->second.x.evaluator->evaluate(1, 0, NULL, NULL, 0, &defaultValue.x);
			values.begin()->second.y.evaluator->evaluate(1, 0, NULL, NULL, 0, &defaultValue.y);
			values.begin()->second.z.evaluator->evaluate(1, 0, NULL, NULL, 0, &defaultValue.z);
			return true;
		}
	}
	return false;
}





