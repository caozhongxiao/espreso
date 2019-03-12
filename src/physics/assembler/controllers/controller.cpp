
#include "esinfo/timeinfo.h"
#include "controller.h"

#include "config/holders/expression.h"

#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/communication.h"
#include "basis/evaluator/evaluator.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;


Controller::Controller(int dimension)
: _dimension(dimension), _dirichletSize(0)
{
	_nDistribution = info::mesh->elements->procNodes->datatarray().distribution();
}

Controller::Parameter::~Parameter()
{
	if (data != NULL) {
		delete data;
	}
}

void Controller::setCoordinates(Parameter &coordinates, serializededata<esint, esint> *procNodes)
{
	if (procNodes == NULL) {
		procNodes = info::mesh->elements->procNodes;
	}

	coordinates.data = new serializededata<esint, double>(coordinates.dimension, procNodes->datatarray().distribution());
	size_t threads = info::env::OMP_NUM_THREADS;

	if (coordinates.dimension == 2) {
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto c = coordinates.data->begin(t);
			for (auto n = procNodes->datatarray().begin(t); n != procNodes->datatarray().end(t); ++n, ++c) {
				c->at(0) = info::mesh->nodes->coordinates->datatarray()[*n].x;
				c->at(1) = info::mesh->nodes->coordinates->datatarray()[*n].y;
			}
		}
	} else {
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto c = coordinates.data->begin(t);
			for (auto n = procNodes->datatarray().begin(t); n != procNodes->datatarray().end(t); ++n, ++c) {
				c->at(0) = info::mesh->nodes->coordinates->datatarray()[*n].x;
				c->at(1) = info::mesh->nodes->coordinates->datatarray()[*n].y;
				c->at(2) = info::mesh->nodes->coordinates->datatarray()[*n].z;
			}
		}
	}
}

void Controller::setDirichlet(std::vector<double> &solution)
{
	std::vector<std::vector<esint> > indices;
	std::vector<double> values;
	dirichletIndices(indices);
	dirichletValues(values);

	for (size_t d = 0, vindex = 0; d < indices.size(); d++) {
		for (size_t i = 0; i < indices[d].size(); i++) {
			solution[indices.size() * indices[d][i] + d] = values[vindex++];
		}
	}
}

void Controller::initKernelParam(Parameter &parameter)
{
	parameter.isConts = false;
	parameter.data = new serializededata<esint, double>(parameter.dimension, info::mesh->elements->procNodes->datatarray().distribution());
}

void Controller::initKernelParam(Parameter &parameter, const std::map<std::string, ECFExpressionVector> &values, double defaultValue)
{
	if (values.size() == 0) { // surely constant
		parameter.isConts = true;
	}
	if (values.size() == 1) { // TODO: generalize for n-dimensinal constant values
	}

	if (parameter.isConts) {
		parameter.data = new serializededata<esint, double>(parameter.dimension, info::mesh->elements->procNodes->datatarray().distribution(), defaultValue);
		return;
	}

	parameter.data = new serializededata<esint, double>(parameter.dimension, info::mesh->elements->procNodes->datatarray().distribution());
}

void Controller::initKernelParam(Parameter &parameter, const std::map<std::string, ECFExpression> &values, double defaultValue)
{
	if (values.size() == 0) { // surely constant
		parameter.isConts = true;
	}
	if (values.size() == 1 && info::mesh->onAllElements(values.begin()->first) && values.begin()->second.evaluator->isConstant()) {
		parameter.isConts = true;
		values.begin()->second.evaluator->evalVector(1, 0, NULL, NULL, 0, &defaultValue);
	}

	if (parameter.isConts) {
		parameter.data = new serializededata<esint, double>(parameter.dimension, info::mesh->elements->procNodes->datatarray().distribution(), defaultValue);
		return;
	}

	parameter.data = new serializededata<esint, double>(parameter.dimension, info::mesh->elements->procNodes->datatarray().distribution());
}

void Controller::initKernelParam(Parameter &parameter, const std::map<std::string, ECFExpression> &values, double defaultValue, BoundaryRegionStore *region)
{
	// TODO: improve (we need to find all ememberships)
	if (values.size() == 0) { // surely constant
		parameter.isConts = true;
	}
	if (values.size() == 1 && info::mesh->onAllElements(values.begin()->first) && values.begin()->second.evaluator->isConstant()) {
		parameter.isConts = true;
		values.begin()->second.evaluator->evalVector(1, 0, NULL, NULL, 0, &defaultValue);
	}

	if (parameter.isConts) {
		parameter.data = new serializededata<esint, double>(parameter.dimension, region->procNodes->datatarray().distribution(), defaultValue);
		return;
	}

	parameter.data = new serializededata<esint, double>(parameter.dimension, region->procNodes->datatarray().distribution());
}

void Controller::initKernelParam(Parameter &parameter, const ECFExpression &value, double defaultValue, BoundaryRegionStore *region)
{
	if (value.evaluator->isConstant()) {
		parameter.isConts = true;
		value.evaluator->evalVector(1, 0, NULL, NULL, 0, &defaultValue);
	}

	if (parameter.isConts) {
		parameter.data = new serializededata<esint, double>(parameter.dimension, region->procNodes->datatarray().distribution(), defaultValue);
		return;
	}

	parameter.data = new serializededata<esint, double>(parameter.dimension, region->procNodes->datatarray().distribution());
}

void Controller::takeKernelParam(Parameter &parameter, Parameter &previous)
{
	parameter = previous;
	parameter.isConts = false; // TODO: check constness
	previous.data = NULL;
}

void Controller::updateKernelParam(Parameter &parameter, const std::map<std::string, ECFExpression> &values, double *cbegin, double *tbegin)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto it = values.begin(); it != values.end(); ++it) {
			ElementsRegionStore *region = info::mesh->eregion(it->first);
			it->second.evaluator->evalFiltered(
					region->elements->datatarray().size(t),
					region->elements->datatarray().begin(t),
					info::mesh->elements->procNodes->boundarytarray().begin(),
					_dimension, cbegin, tbegin, time::current,
					parameter.data->datatarray().data()
			);
		}
	}
}

void Controller::updateKernelParam(Parameter &parameter, const std::map<std::string, ECFExpressionVector> &values, double *cbegin, double *tbegin)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto it = values.begin(); it != values.end(); ++it) {
			ElementsRegionStore *region = info::mesh->eregion(it->first);
			it->second.x.evaluator->evalFiltered(
					region->elements->datatarray().size(t), _dimension,
					region->elements->datatarray().begin(t),
					info::mesh->elements->procNodes->boundarytarray().begin(),
					_dimension, cbegin, tbegin, time::current,
					parameter.data->datatarray().data() + 0
			);
			it->second.y.evaluator->evalFiltered(
					region->elements->datatarray().size(t), _dimension,
					region->elements->datatarray().begin(t),
					info::mesh->elements->procNodes->boundarytarray().begin(),
					_dimension, cbegin, tbegin, time::current,
					parameter.data->datatarray().data() + 1
			);
			if (_dimension == 3) {
				it->second.z.evaluator->evalFiltered(
						region->elements->datatarray().size(t), _dimension,
						region->elements->datatarray().begin(t),
						info::mesh->elements->procNodes->boundarytarray().begin(),
						_dimension, cbegin, tbegin, time::current,
						parameter.data->datatarray().data() + 2
				);
			}
		}
	}
}

void Controller::updateKernelParam(Parameter &parameter, const ECFExpression &value, double *cbegin, double *tbegin,
		BoundaryRegionStore *region)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t offset = parameter.data->datatarray().distribution()[t];
		double *_cbegin = cbegin != NULL ? cbegin + _dimension * offset : NULL;
		double *_tbegin = tbegin != NULL ? tbegin + offset : NULL;
		value.evaluator->evalVector(
				parameter.data->datatarray().size(t),
				_dimension, _cbegin, _tbegin, time::current,
				parameter.data->datatarray().data() + parameter.dimension * offset
		);
	}
}

void Controller::processBEMdomain(esint domain, double *values)
{
	eslog::error("ESPRESO internal error: BEM is not implemented for the requested physics.\n");
}

void Controller::fillBEMInterior(esint domain, double *values)
{
	eslog::error("ESPRESO internal error: BEM is not implemented for the requested physics.\n");
}

void Controller::kernelToBoundary(Parameter &parameter, Parameter &boundary, BoundaryRegionStore *region)
{
	eslog::error("ESPRESO internal error: kernelToBoundary.\n");
}

void Controller::kernelToNodes(Parameter &kvalues, std::vector<double> &nvalues)
{
	std::fill(nvalues.begin(), nvalues.end(), 0);
	auto i = kvalues.data->datatarray().begin();
	for (auto n = info::mesh->elements->procNodes->datatarray().cbegin(); n != info::mesh->elements->procNodes->datatarray().cend(); ++n, ++i) {
		for (int d = 0; d < kvalues.dimension; d++) {
			nvalues[*n * kvalues.dimension + d] += *i;
		}
	}

	auto &nelements = info::mesh->nodes->elements->boundarytarray();
	for (size_t i = 0; i < nvalues.size() / kvalues.dimension; i++) {
		for (int d = 0; d < kvalues.dimension; d++) {
			nvalues[i * kvalues.dimension + d] /= nelements[i + 1] - nelements[i];
		}
	}

	std::vector<std::vector<double> > sBuffer(info::mesh->neighbours.size()), rBuffer(info::mesh->neighbours.size());

	auto nranks = info::mesh->nodes->ranks->begin();
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		if (nranks->size() > 1) {
			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != info::mpi::rank) {
					while (info::mesh->neighbours[noffset] < *r) {
						++noffset;
					}
					for (int d = 0; d < kvalues.dimension; d++) {
						sBuffer[noffset].push_back(nvalues[n * kvalues.dimension + d]);
					}
				}
			}
		}
	}

	for (size_t n = 0; n < info::mesh->neighbours.size(); n++) {
		rBuffer[n].resize(sBuffer[n].size());
	}

	if (!Communication::exchangeKnownSize(sBuffer, rBuffer, info::mesh->neighbours)) {
		eslog::error("ESPRESO internal error: exchange diagonal values.\n");
	}

	nranks = info::mesh->nodes->ranks->begin();
	std::vector<esint> nindex(info::mesh->neighbours.size());
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		if (nranks->size() > 1) {
			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != info::mpi::rank) {
					while (info::mesh->neighbours[noffset] < *r) {
						++noffset;
					}
					for (int d = 0; d < kvalues.dimension; d++) {
						nvalues[n * kvalues.dimension + d] += rBuffer[noffset][nindex[noffset]++];
					}
				}
			}
		}
	}
}

void Controller::kernelToElements(Parameter &kvalues, std::vector<double> &evalues)
{
	size_t threads = info::env::OMP_NUM_THREADS;
	auto elements = info::mesh->elements;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t noffset = elements->procNodes->datatarray().distribution()[t];
		size_t eoffset = elements->distribution[t];
		for (auto enodes = elements->procNodes->cbegin(t); enodes != elements->procNodes->cend(t); ++enodes, ++eoffset) {
			for (int d = 0; d < kvalues.dimension; d++) {
				double sum = 0;
				for (auto n = enodes->begin(); n != enodes->end(); ++n, ++noffset) {
					sum += kvalues.data->datatarray()[kvalues.dimension * noffset + d];
				}
				evalues[kvalues.dimension * eoffset + d] = sum / enodes->size();
				noffset -= enodes->size();
			}
			noffset += enodes->size();
		}
	}
}

void Controller::nodesToKernels(std::vector<double> &nvalues, Parameter &kvalues, serializededata<esint, esint> *procNodes)
{
	size_t threads = info::env::OMP_NUM_THREADS;
	if (procNodes == NULL) {
		procNodes = info::mesh->elements->procNodes;
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto i = kvalues.data->datatarray().begin(t);
		for (auto n = procNodes->datatarray().cbegin(t); n != procNodes->datatarray().cend(t); ++n) {
			for (int d = 0; d < kvalues.dimension; ++d, ++i) {
				*i = nvalues[kvalues.dimension * *n + d];
			}
		}
	}
}





