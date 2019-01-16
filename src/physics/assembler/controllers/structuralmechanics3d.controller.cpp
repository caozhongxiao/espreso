
#include "physics/assembler/dataholder.h"
#include "esinfo/time.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "structuralmechanics3d.controller.h"
#include "physics/assembler/kernels/structuralmechanics3d.kernel.h"

#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "config/ecf/root.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"

using namespace espreso;

StructuralMechanics3DControler::StructuralMechanics3DControler(StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanicsControler(configuration)
{
	_kernel = new StructuralMechanics3DKernel();

	double defaultTemperature = 0;

	_ncoordinate.data = new serializededata<esint, double>(3, _nDistribution);
	_ntemperature.data = new serializededata<esint, double>(1, _nDistribution);

	_nInitialTemperature.isConts = setDefault(info::ecf->structural_mechanics_3d.initial_temperature, defaultTemperature);
	_nInitialTemperature.data = new serializededata<esint, double>(1, _nDistribution, defaultTemperature);

	_nacceleration.isConts = false;
	_nacceleration.data = new serializededata<esint, double>(3, _nDistribution);

	_nangularVelocity.isConts = false;
	_nangularVelocity.data = new serializededata<esint, double>(3, _nDistribution);

	_displacement = info::mesh->nodes->appendData(3, { "DISPLACEMENT", "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z" });

	_boundaries.resize(info::mesh->boundaryRegions.size());
}

StructuralMechanics3DControler::~StructuralMechanics3DControler()
{
	delete _kernel;
}

void StructuralMechanics3DControler::dirichletIndices(std::vector<std::vector<esint> > &indices)
{
	indices.resize(3);

	for (auto it = _configuration.displacement.regions.begin(); it != _configuration.displacement.regions.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		if (it->second.all.value.size() || it->second.x.value.size()) {
			indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
		}
		if (it->second.all.value.size() || it->second.y.value.size()) {
			indices[1].insert(indices[1].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
		}
		if (it->second.all.value.size() || it->second.z.value.size()) {
			indices[2].insert(indices[2].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
		}
	}

	for (auto it = _configuration.displacement.intersections.begin(); it != _configuration.displacement.intersections.end(); ++it) {
		BoundaryRegionsIntersectionStore *region = info::mesh->ibregion(it->first);
		if (it->second.all.value.size() || it->second.x.value.size()) {
			indices[0].insert(indices[0].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
		}
		if (it->second.all.value.size() || it->second.y.value.size()) {
			indices[1].insert(indices[1].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
		}
		if (it->second.all.value.size() || it->second.z.value.size()) {
			indices[2].insert(indices[2].end(), region->uniqueNodes->datatarray().begin(), region->uniqueNodes->datatarray().end());
		}
	}
	_dirichletSize = indices[0].size() + indices[1].size() + indices[2].size();
}

void StructuralMechanics3DControler::dirichletValues(std::vector<double> &values)
{
	values.resize(_dirichletSize);

	size_t offset = 0;
	double *coors = reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data());
	auto eval = [&] (Evaluator *evaluator, tarray<esint> &nodes) {
		evaluator->evalSelected(nodes.size(), nodes.data(), 3, coors, NULL, time::current, values.data() + offset);
		offset += nodes.size();
	};

	auto pick = [&] (ECFExpressionOptionalVector &vector, tarray<esint> &nodes) {
		if (vector.all.value.size()) {
			eval(vector.all.evaluator, nodes);
			eval(vector.all.evaluator, nodes);
			eval(vector.all.evaluator, nodes);
		} else {
			if (vector.x.value.size()) {
				eval(vector.x.evaluator, nodes);
			}
			if (vector.y.value.size()) {
				eval(vector.y.evaluator, nodes);
			}
			if (vector.z.value.size()) {
				eval(vector.z.evaluator, nodes);
			}
		}
	};

	for (auto it = _configuration.displacement.regions.begin(); it != _configuration.displacement.regions.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		pick(it->second, region->uniqueNodes->datatarray());
	}

	for (auto it = _configuration.displacement.intersections.begin(); it != _configuration.displacement.intersections.end(); ++it) {
		BoundaryRegionsIntersectionStore *region = info::mesh->ibregion(it->first);
		pick(it->second, region->uniqueNodes->datatarray());
	}
}

void StructuralMechanics3DControler::initData()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto c = _ncoordinate.data->begin(t);
		for (auto n = info::mesh->elements->procNodes->datatarray().begin(t); n != info::mesh->elements->procNodes->datatarray().end(t); ++n, ++c) {
			c->at(0) = info::mesh->nodes->coordinates->datatarray()[*n].x;
			c->at(1) = info::mesh->nodes->coordinates->datatarray()[*n].y;
			c->at(2) = info::mesh->nodes->coordinates->datatarray()[*n].z;
		}
	}

	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = time::current;

	updateERegions(info::ecf->structural_mechanics_3d.initial_temperature, _nInitialTemperature.data->datatarray(), 1, cbegin, tbegin, time);
	updateERegions(_configuration.acceleration, _nacceleration.data->datatarray(), 3, cbegin, tbegin, time);
	updateERegions(_configuration.angular_velocity, _nangularVelocity.data->datatarray(), 3, cbegin, tbegin, time);

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 2) { // TODO: implement edge processing

			auto &distribution = region->procNodes->datatarray().distribution();

			_boundaries[r].coordinate.data = new serializededata<esint, double>(3, distribution);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto c = _boundaries[r].coordinate.data->begin(t);
				for (auto n = region->procNodes->datatarray().begin(t); n != region->procNodes->datatarray().end(t); ++n, ++c) {
					c->at(0) = info::mesh->nodes->coordinates->datatarray()[*n].x;
					c->at(1) = info::mesh->nodes->coordinates->datatarray()[*n].y;
					c->at(2) = info::mesh->nodes->coordinates->datatarray()[*n].z;
				}
			}

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();

			auto pressure = _configuration.normal_pressure.find(region->name);
			if (pressure != _configuration.normal_pressure.end()) {
				updateBRegions(pressure->second, _boundaries[r].normalPressure, distribution, 3, cbegin, tbegin, time);
			}
		}
	}
}

void StructuralMechanics3DControler::nextTime()
{
	if (time::isInitial()) {
		return;
	}

	parametersChanged();
}

void StructuralMechanics3DControler::parametersChanged()
{
	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = time::current;

	updateERegions(_configuration.acceleration, _nacceleration.data->datatarray(), 3, cbegin, tbegin, time);
	updateERegions(_configuration.angular_velocity, _nangularVelocity.data->datatarray(), 3, cbegin, tbegin, time);

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 2) {

			auto &distribution = region->procNodes->datatarray().distribution();

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();

			auto pressure = _configuration.normal_pressure.find(region->name);
			if (pressure != _configuration.normal_pressure.end()) {
				updateBRegions(pressure->second, _boundaries[r].normalPressure, distribution, 2, cbegin, tbegin, time);
			}
		}
	}
}

void StructuralMechanics3DControler::processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler)
{
	auto enodes = info::mesh->elements->procNodes->cbegin() + filler.begin;
	StructuralMechanics3DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin();
	iterator.temperature        = _ntemperature.data->datatarray().begin() + noffset;
	iterator.initialTemperature = _nInitialTemperature.data->datatarray().begin() + noffset;
	iterator.coordinates        = _ncoordinate.data->datatarray().begin() + noffset * 3;
	iterator.acceleration       = _nacceleration.data->datatarray().begin() + noffset * 3;
	iterator.angularVelocity    = _nangularVelocity.data->datatarray().begin() + noffset * 3;


	for (esint e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = info::mesh->elements->epointers->datatarray()[e];
		iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

		_kernel->processElement(matrices, parameters, iterator, filler.Ke, filler.Me, filler.Re, filler.fe);
		filler.insert(3 * enodes->size());

		iterator.temperature        += enodes->size();
		iterator.initialTemperature += enodes->size();
		iterator.coordinates        += enodes->size() * 3;
		iterator.acceleration       += enodes->size() * 3;
		iterator.angularVelocity    += enodes->size() * 3;
	}
}

void StructuralMechanics3DControler::processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler)
{
	if (info::mesh->boundaryRegions[rindex]->dimension != 2) {
		return;
	}

	auto enodes = info::mesh->boundaryRegions[rindex]->procNodes->cbegin() + filler.begin;
	StructuralMechanics3DKernel::BoundaryIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->boundaryRegions[rindex]->procNodes->datatarray().begin();
	iterator.coordinates = _boundaries[rindex].coordinate.data->datatarray().begin() + noffset * 3;

	iterator.normalPressure = _boundaries[rindex].normalPressure.data ? _boundaries[rindex].normalPressure.data->datatarray().begin() + noffset : NULL;

	for (esint e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = info::mesh->boundaryRegions[rindex]->epointers->datatarray()[e];

		_kernel->processFace(matrices, parameters, iterator, filler.Ke, filler.fe);
		filler.insert(3 * enodes->size());

		if (iterator.normalPressure) {
			iterator.normalPressure += enodes->size();
		}
	}
}

void StructuralMechanics3DControler::processSolution()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		auto enodes = info::mesh->elements->procNodes->cbegin(t);
		StructuralMechanics3DKernel::SolutionIterator iterator;

		size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin(t);
		iterator.temperature = _ntemperature.data->datatarray().begin(t) + noffset;
		iterator.coordinates = _ncoordinate.data->datatarray().begin(t) + noffset * 3;

		for (size_t e = info::mesh->elements->distribution[t]; e < info::mesh->elements->distribution[t + 1]; ++e, ++enodes) {
			iterator.element = info::mesh->elements->epointers->datatarray()[e];
			iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

			_kernel->processSolution(iterator);

			iterator.temperature += enodes->size();
			iterator.coordinates += enodes->size() * 3;
		}
	}
}







