
#include "controller.h"

#include "../../../config/ecf/environment.h"

#include "../../../globals/run.h"
#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../config/expression.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/elementsregionstore.h"

using namespace espreso;


Controler::Controler()
: _dirichletSize(0)
{
	_nDistribution = run::mesh->elements->procNodes->datatarray().distribution();
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
	if (run::mesh->onAllElements(values.begin()->first) && values.begin()->second.evaluator->isConstant()) {
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
	if (run::mesh->onAllElements(values.begin()->first)) {
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
		eslocal csize, double *cbegin, double *tbegin, double time,
		EvaluatorParameters updatedParams)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto it = settings.begin(); it != settings.end(); ++it) {
			ElementsRegionStore *region = run::mesh->eregion(it->first);
			it->second.evaluator->evalFiltered(
					region->elements->datatarray().size(t),
					region->elements->datatarray().begin(t),
					run::mesh->elements->procNodes->boundarytarray().begin(),
					csize, cbegin, tbegin, time, data.data()
			);
		}
	}
}

void Controler::updateERegions(
		const std::map<std::string, ECFExpressionVector> &settings, tarray<double> &data,
		eslocal csize, double *cbegin, double *tbegin, double time,
		EvaluatorParameters updatedParams)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto it = settings.begin(); it != settings.end(); ++it) {
			ElementsRegionStore *region = run::mesh->eregion(it->first);
			it->second.x.evaluator->evalFiltered(
					region->elements->datatarray().size(t), csize,
					region->elements->datatarray().begin(t),
					run::mesh->elements->procNodes->boundarytarray().begin(t),
					csize, cbegin, tbegin, time, data.data() + 0
			);
			it->second.y.evaluator->evalFiltered(
					region->elements->datatarray().size(t), csize,
					region->elements->datatarray().begin(t),
					run::mesh->elements->procNodes->boundarytarray().begin(t),
					csize, cbegin, tbegin, time, data.data() + 1
			);
			if (csize == 3) {
				it->second.z.evaluator->evalFiltered(
						region->elements->datatarray().size(t), csize,
						region->elements->datatarray().begin(t),
						run::mesh->elements->procNodes->boundarytarray().begin(t),
						csize, cbegin, tbegin, time, data.data() + 2
				);
			}
		}
	}
}

void Controler::updateBRegions(
		const ECFExpression &expression, Parameter &parameter, const std::vector<size_t> &distribution,
		eslocal csize, double *cbegin, double *tbegin, double time,
		EvaluatorParameters updatedParams)
{
	size_t threads = environment->OMP_NUM_THREADS;

	parameter.data = new serializededata<eslocal, double>(1, distribution);
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
	for (auto n = run::mesh->elements->procNodes->datatarray().cbegin(); n != run::mesh->elements->procNodes->datatarray().cend(); ++n, ++i) {
		averagedData[*n] += *i;
	}

	auto &nelements = run::mesh->nodes->elements->boundarytarray();
	for (size_t i = 0; i < averagedData.size(); i++) {
		averagedData[i] /= nelements[i + 1] - nelements[i];
	}
}

void Controler::nodeValuesToElements(tarray<double> &nodeData, std::vector<double> &elementData)
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t noffset = run::mesh->elements->procNodes->cbegin(t)->begin() - run::mesh->elements->procNodes->cbegin()->begin();
		size_t eoffset = run::mesh->elements->distribution[t];
		for (auto enodes = run::mesh->elements->procNodes->cbegin(t); enodes != run::mesh->elements->procNodes->cend(t); ++enodes) {
			double sum = 0;
			for (auto n = enodes->begin(); n != enodes->end(); ++n, ++noffset) {
				sum += nodeData[noffset];
			}
			elementData[eoffset] = sum / enodes->size();
		}
	}
}





