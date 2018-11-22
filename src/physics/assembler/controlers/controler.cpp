
#include "controler.h"

#include "../../../config/ecf/environment.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../config/expression.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/elementsregionstore.h"

using namespace espreso;


Controler::Controler(Mesh &mesh, const Step &step)
: _mesh(mesh), _step(step), _dirichletSize(0)
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
		values.begin()->second.evaluator->evalVector(1, 0, NULL, NULL, 0, &defaultValue);
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

			values.begin()->second.x.evaluator->evalVector(1, 0, NULL, NULL, 0, &defaultValue.x);
			values.begin()->second.y.evaluator->evalVector(1, 0, NULL, NULL, 0, &defaultValue.y);
			values.begin()->second.z.evaluator->evalVector(1, 0, NULL, NULL, 0, &defaultValue.z);
			return true;
		}
	}
	return false;
}

void Controler::evaluate(
		const std::map<std::string, ECFExpression> &settings, tarray<double> &data,
		eslocal csize, double *cbegin)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto it = settings.begin(); it != settings.end(); ++it) {
			ElementsRegionStore *region = _mesh.eregion(it->first);
			it->second.evaluator->evalFiltered(
					region->elements->datatarray().size(t),
					region->elements->datatarray().begin(t),
					_mesh.elements->procNodes->boundarytarray().begin(t),
					2, cbegin, NULL, 0, data.data()
			);
		}
	}
}

void Controler::evaluate(
		const std::map<std::string, ECFExpressionVector> &settings, tarray<double> &data,
		eslocal csize, double *cbegin)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto it = settings.begin(); it != settings.end(); ++it) {
			ElementsRegionStore *region = _mesh.eregion(it->first);
			it->second.x.evaluator->evalFiltered(
					region->elements->datatarray().size(t), 2,
					region->elements->datatarray().begin(t),
					_mesh.elements->procNodes->boundarytarray().begin(t),
					2, cbegin, NULL, 0, data.data()
			);
			it->second.y.evaluator->evalFiltered(
					region->elements->datatarray().size(t), 2,
					region->elements->datatarray().begin(t),
					_mesh.elements->procNodes->boundarytarray().begin(t),
					2, cbegin, NULL, 0, data.data() + 1
			);
		}
	}
}

void Controler::nodeValuesToElements(tarray<double> &nodeData, std::vector<double> &elementData)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t noffset = _mesh.elements->procNodes->cbegin(t)->begin() - _mesh.elements->procNodes->cbegin()->begin();
		size_t eoffset = _mesh.elements->distribution[t];
		for (auto enodes = _mesh.elements->procNodes->cbegin(t); enodes != _mesh.elements->procNodes->cend(t); ++enodes) {
			double sum = 0;
			for (auto n = enodes->begin(); n != enodes->end(); ++n, ++noffset) {
				sum += nodeData[noffset];
			}
			elementData[eoffset] = sum / enodes->size();
		}
	}
}





