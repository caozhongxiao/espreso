
#include "config/holders/expression.h"
#include "esinfo/meshinfo.h"
#include "controller.h"

#include "config/ecf/environment.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"

using namespace espreso;


Controler::Controler()
: _dirichletSize(0)
{
	_nDistribution = info::mesh->elements->procNodes->datatarray().distribution();
}

Controler::Parameter::~Parameter()
{
	if (data != NULL) {
		delete data;
	}
}

bool Controler::setDefault(const std::map<std::string, ECFExpression> &values, double &defaultValue)
{
	if (values.size() == 0) {
		return true;
	}
	if (values.size() > 1) {
		return false;
	}
	if (info::mesh->onAllElements(values.begin()->first) && values.begin()->second.evaluator->isConstant()) {
		values.begin()->second.evaluator->evalVector(1, 0, NULL, NULL, 0, &defaultValue);
		return true;
	}
	return false;
}

bool Controler::setDefault(const std::map<std::string, ECFExpressionVector> &values, Point &defaultValue)
{
	if (values.size() == 0) {
		return true;
	}
	if (values.size() > 1) {
		return false;
	}
	if (info::mesh->onAllElements(values.begin()->first)) {
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

void Controler::updateERegions(
		const std::map<std::string, ECFExpression> &settings, tarray<double> &data,
		esint csize, double *cbegin, double *tbegin, double time,
		EvaluatorParameters updatedParams)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto it = settings.begin(); it != settings.end(); ++it) {
			ElementsRegionStore *region = info::mesh->eregion(it->first);
			it->second.evaluator->evalFiltered(
					region->elements->datatarray().size(t),
					region->elements->datatarray().begin(t),
					info::mesh->elements->procNodes->boundarytarray().begin(),
					csize, cbegin, tbegin, time, data.data()
			);
		}
	}
}

void Controler::updateERegions(
		const std::map<std::string, ECFExpressionVector> &settings, tarray<double> &data,
		esint csize, double *cbegin, double *tbegin, double time,
		EvaluatorParameters updatedParams)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto it = settings.begin(); it != settings.end(); ++it) {
			ElementsRegionStore *region = info::mesh->eregion(it->first);
			it->second.x.evaluator->evalFiltered(
					region->elements->datatarray().size(t), csize,
					region->elements->datatarray().begin(t),
					info::mesh->elements->procNodes->boundarytarray().begin(),
					csize, cbegin, tbegin, time, data.data() + 0
			);
			it->second.y.evaluator->evalFiltered(
					region->elements->datatarray().size(t), csize,
					region->elements->datatarray().begin(t),
					info::mesh->elements->procNodes->boundarytarray().begin(),
					csize, cbegin, tbegin, time, data.data() + 1
			);
			if (csize == 3) {
				it->second.z.evaluator->evalFiltered(
						region->elements->datatarray().size(t), csize,
						region->elements->datatarray().begin(t),
						info::mesh->elements->procNodes->boundarytarray().begin(),
						csize, cbegin, tbegin, time, data.data() + 2
				);
			}
		}
	}
}

void Controler::updateBRegions(
		const ECFExpression &expression, Parameter &parameter, const std::vector<size_t> &distribution,
		esint csize, double *cbegin, double *tbegin, double time,
		EvaluatorParameters updatedParams)
{
	size_t threads = environment->OMP_NUM_THREADS;

	parameter.data = new serializededata<esint, double>(1, distribution);
	parameter.isConts = expression.evaluator->isConstant();

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		cbegin += distribution[t] * csize;
		tbegin += distribution[t];
		double *tresults = parameter.data->datatarray().begin(t);

		expression.evaluator->evalVector(parameter.data->datatarray().size(t), 2, cbegin, tbegin, time, tresults);
	}
}

void Controler::averageNodeInitilization(tarray<double> &initData, std::vector<double> &averagedData)
{
	auto i = initData.begin();
	for (auto n = info::mesh->elements->procNodes->datatarray().cbegin(); n != info::mesh->elements->procNodes->datatarray().cend(); ++n, ++i) {
		averagedData[*n] += *i;
	}

	auto &nelements = info::mesh->nodes->elements->boundarytarray();
	for (size_t i = 0; i < averagedData.size(); i++) {
		averagedData[i] /= nelements[i + 1] - nelements[i];
	}
}

void Controler::nodeValuesToElements(int dimension, tarray<double> &nodeData, std::vector<double> &elementData)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t noffset = info::mesh->elements->procNodes->cbegin(t)->begin() - info::mesh->elements->procNodes->cbegin()->begin();
		size_t eoffset = info::mesh->elements->distribution[t];
		for (auto enodes = info::mesh->elements->procNodes->cbegin(t); enodes != info::mesh->elements->procNodes->cend(t); ++enodes, ++eoffset) {
			for (int d = 0; d < dimension; d++) {
				double sum = 0;
				for (auto n = enodes->begin(); n != enodes->end(); ++n, ++noffset) {
					sum += nodeData[dimension * noffset + d];
				}
				elementData[dimension * eoffset + d] = sum / enodes->size();
				noffset -= enodes->size();
			}
			noffset += enodes->size();
		}
	}
}





