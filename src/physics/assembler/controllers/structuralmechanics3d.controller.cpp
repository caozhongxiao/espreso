
#include "esinfo/time.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "structuralmechanics3d.controller.h"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/kernels/structuralmechanics3d.kernel.h"

#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "config/ecf/root.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

StructuralMechanics3DController::StructuralMechanics3DController(StructuralMechanics3DController *previous, StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanicsController(3, previous, configuration), _kernel(new StructuralMechanics3DKernel())
{
	if (previous) {
		takeKernelParam(_kcoordinate, previous->_kcoordinate);
		takeKernelParam(_kinitialTemperature, previous->_kinitialTemperature);
	} else {
		_ndisplacement = info::mesh->nodes->appendData(3, { "DISPLACEMENT", "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z" });

		setCoordinates(_kcoordinate);
		initKernelParam(_kinitialTemperature, info::ecf->structural_mechanics_3d.initial_temperature, 0);
		initKernelParam(_kthickness, info::ecf->structural_mechanics_3d.thickness, 1);

		if (!_kinitialTemperature.isConts) {
			updateKernelParam(_kinitialTemperature, info::ecf->structural_mechanics_3d.initial_temperature, _kcoordinate.data->datatarray().data(), NULL);
		}
	}

	initKernelParam(_ktemperature, configuration.temperature, 0);
	initKernelParam(_kacceleration, configuration.acceleration, 0);
	initKernelParam(_kangularVelocity, configuration.angular_velocity, 0);

	_boundaries.resize(info::mesh->boundaryRegions.size(), BoundaryParameters(3));
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 2) {
			setCoordinates(_boundaries[r].coordinate, region->procNodes);

			auto pressure = _configuration.normal_pressure.find(region->name);
			if (pressure != _configuration.normal_pressure.end()) {
				initKernelParam(_boundaries[r].normalPressure, pressure->second, 0, region);
			}
		}
	}
}

StructuralMechanics3DController::~StructuralMechanics3DController()
{
	if (_kernel) {
		delete _kernel;
	}
}

const PhysicsConfiguration& StructuralMechanics3DController::configuration() const
{
	return info::ecf->structural_mechanics_3d;
}

void StructuralMechanics3DController::dirichletIndices(std::vector<std::vector<esint> > &indices)
{
	indices.resize(3);

	for (auto it = _configuration.displacement.begin(); it != _configuration.displacement.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		if (it->second.all.value.size() || it->second.x.value.size()) {
			indices[0].insert(indices[0].end(), region->nodes->datatarray().begin(), region->nodes->datatarray().end());
		}
		if (it->second.all.value.size() || it->second.y.value.size()) {
			indices[1].insert(indices[1].end(), region->nodes->datatarray().begin(), region->nodes->datatarray().end());
		}
		if (it->second.all.value.size() || it->second.z.value.size()) {
			indices[2].insert(indices[2].end(), region->nodes->datatarray().begin(), region->nodes->datatarray().end());
		}
	}
	_dirichletSize = indices[0].size() + indices[1].size() + indices[2].size();
}

void StructuralMechanics3DController::dirichletValues(std::vector<double> &values)
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

	for (auto it = _configuration.displacement.begin(); it != _configuration.displacement.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		pick(it->second, region->nodes->datatarray());
	}
}

void StructuralMechanics3DController::nextTime()
{
	parametersChanged();
}

void StructuralMechanics3DController::parametersChanged()
{
	double *cbegin = _kcoordinate.data->datatarray().data();
	double *tbegin = _ktemperature.data->datatarray().data();

	if (!_ktemperature.isConts) {
		updateKernelParam(_ktemperature, _configuration.temperature, cbegin, NULL);
	}
	if (!_kacceleration.isConts) {
		updateKernelParam(_kacceleration, _configuration.acceleration, cbegin, tbegin);
	}
	if (!_kangularVelocity.isConts) {
		updateKernelParam(_kangularVelocity, _configuration.angular_velocity, cbegin, tbegin);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 2) {
			cbegin = _boundaries[r].coordinate.data->datatarray().begin();

			auto pressure = _configuration.normal_pressure.find(region->name);
			if (pressure != _configuration.normal_pressure.end()) {
				if (!_boundaries[r].normalPressure.isConts) {
					updateKernelParam(_boundaries[r].normalPressure, pressure->second, cbegin, NULL, region);
				}
			}
		}
	}
}

void StructuralMechanics3DController::processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler)
{
	auto enodes = info::mesh->elements->procNodes->cbegin() + filler.begin;
	StructuralMechanics3DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin();
	iterator.temperature        = _ktemperature.data->datatarray().begin() + noffset;
	iterator.initialTemperature = _kinitialTemperature.data->datatarray().begin() + noffset;
	iterator.coordinates        = _kcoordinate.data->datatarray().begin() + noffset * 3;
	iterator.acceleration       = _kacceleration.data->datatarray().begin() + noffset * 3;
	iterator.angularVelocity    = _kangularVelocity.data->datatarray().begin() + noffset * 3;


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

void StructuralMechanics3DController::processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler)
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

void StructuralMechanics3DController::processSolution()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		auto enodes = info::mesh->elements->procNodes->cbegin(t);
		StructuralMechanics3DKernel::SolutionIterator iterator;

		size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin(t);
		iterator.temperature = _ktemperature.data->datatarray().begin(t) + noffset;
		iterator.coordinates = _kcoordinate.data->datatarray().begin(t) + noffset * 3;

		for (size_t e = info::mesh->elements->distribution[t]; e < info::mesh->elements->distribution[t + 1]; ++e, ++enodes) {
			iterator.element = info::mesh->elements->epointers->datatarray()[e];
			iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

			_kernel->processSolution(iterator);

			iterator.temperature += enodes->size();
			iterator.coordinates += enodes->size() * 3;
		}
	}
}







