
#include "esinfo/time.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "structuralmechanics2d.controller.h"
#include "physics/assembler/dataholder.h"
#include "physics/assembler/kernels/structuralmechanics2d.kernel.h"

#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "config/ecf/root.h"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

StructuralMechanics2DController::StructuralMechanics2DController(StructuralMechanics2DController *previous, StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanicsController(2, previous, configuration), _kernel(new StructuralMechanics2DKernel())
{
	if (previous) {
		takeKernelParam(_kcoordinate, previous->_kcoordinate);
		takeKernelParam(_kinitialTemperature, previous->_kinitialTemperature);
		takeKernelParam(_kthickness, previous->_kthickness);
	} else {
		_ndisplacement = info::mesh->nodes->appendData(2, { "DISPLACEMENT", "DISPLACEMENT_X", "DISPLACEMENT_Y" });
		if (info::ecf->output.results_selection.thickness) {
			_ethickness = info::mesh->elements->appendData(1, { "THICKNESS" });
		}

		setCoordinates(_kcoordinate);
		initKernelParam(_kinitialTemperature, info::ecf->structural_mechanics_2d.initial_temperature, 0);
		initKernelParam(_kthickness, info::ecf->structural_mechanics_2d.thickness, 1);

		if (!_kinitialTemperature.isConts) {
			updateKernelParam(_kinitialTemperature, info::ecf->structural_mechanics_2d.initial_temperature, _kcoordinate.data->datatarray().data(), NULL);
		}
		if (_ethickness && _kthickness.isConts) {
			kernelToElements(_kthickness, _ethickness->data);
		}
	}

	initKernelParam(_ktemperature, configuration.temperature, 0);
	initKernelParam(_kacceleration, configuration.acceleration, 0);
	initKernelParam(_kangularVelocity, configuration.angular_velocity, 0);

	_boundaries.resize(info::mesh->boundaryRegions.size(), BoundaryParameters(2));
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 1) {
			setCoordinates(_boundaries[r].coordinate, region->procNodes);
			initKernelParam(_boundaries[r].thickness, info::ecf->structural_mechanics_2d.thickness, 1, region);

			auto pressure = _configuration.normal_pressure.find(region->name);
			if (pressure != _configuration.normal_pressure.end()) {
				initKernelParam(_boundaries[r].normalPressure, pressure->second, 0, region);
			}
		}
	}
}

StructuralMechanics2DController::~StructuralMechanics2DController()
{
	if (_kernel) {
		delete _kernel;
	}
}

const PhysicsConfiguration& StructuralMechanics2DController::configuration() const
{
	return info::ecf->structural_mechanics_2d;
}


void StructuralMechanics2DController::dirichletIndices(std::vector<std::vector<esint> > &indices)
{
	indices.resize(2);

	for (auto it = _configuration.displacement.begin(); it != _configuration.displacement.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		if (it->second.all.value.size() || it->second.x.value.size()) {
			indices[0].insert(indices[0].end(), region->nodes->datatarray().begin(), region->nodes->datatarray().end());
		}
		if (it->second.all.value.size() || it->second.y.value.size()) {
			indices[1].insert(indices[1].end(), region->nodes->datatarray().begin(), region->nodes->datatarray().end());
		}
	}
	_dirichletSize = indices[0].size() + indices[1].size();
}

void StructuralMechanics2DController::dirichletValues(std::vector<double> &values)
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
		} else {
			if (vector.x.value.size()) {
				eval(vector.x.evaluator, nodes);
			}
			if (vector.y.value.size()) {
				eval(vector.y.evaluator, nodes);
			}
		}
	};

	for (auto it = _configuration.displacement.begin(); it != _configuration.displacement.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		pick(it->second, region->nodes->datatarray());
	}
}

void StructuralMechanics2DController::nextTime()
{
	parametersChanged();
}

void StructuralMechanics2DController::parametersChanged()
{
	double *cbegin = _kcoordinate.data->datatarray().data();
	double *tbegin = _ktemperature.data->datatarray().data();

	if (!_ktemperature.isConts) {
		updateKernelParam(_ktemperature, _configuration.temperature, cbegin, NULL);
	}
	if (!_kthickness.isConts) {
		updateKernelParam(_kthickness, info::ecf->structural_mechanics_2d.thickness, cbegin, tbegin);
		if (_ethickness) {
			kernelToElements(_kthickness, _ethickness->data);
		}
	}
	if (!_kacceleration.isConts) {
		updateKernelParam(_kacceleration, _configuration.acceleration, cbegin, tbegin);
	}
	if (!_kangularVelocity.isConts) {
		updateKernelParam(_kangularVelocity, _configuration.angular_velocity, cbegin, tbegin);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 1) {
			cbegin = _boundaries[r].coordinate.data->datatarray().begin();
			if (!_boundaries[r].thickness.isConts) {
				kernelToBoundary(_kthickness, _boundaries[r].thickness, region);
			}

			auto pressure = _configuration.normal_pressure.find(region->name);
			if (pressure != _configuration.normal_pressure.end()) {
				if (!_boundaries[r].normalPressure.isConts) {
					updateKernelParam(_boundaries[r].normalPressure, pressure->second, cbegin, NULL, region);
				}
			}
		}
	}
}

void StructuralMechanics2DController::processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler)
{
	auto enodes = info::mesh->elements->procNodes->cbegin() + filler.begin;
	StructuralMechanics2DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin();
	iterator.temperature        = _ktemperature.data->datatarray().begin() + noffset;
	iterator.initialTemperature = _kinitialTemperature.data->datatarray().begin() + noffset;
	iterator.coordinates        = _kcoordinate.data->datatarray().begin() + noffset * 2;
	iterator.acceleration       = _kacceleration.data->datatarray().begin() + noffset * 2;
	iterator.angularVelocity    = _kangularVelocity.data->datatarray().begin() + noffset * 3;
	iterator.thickness          = _kthickness.data->datatarray().begin() + noffset;

	for (esint e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = info::mesh->elements->epointers->datatarray()[e];
		iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

		_kernel->processElement(matrices, parameters, iterator, filler.Ke, filler.Me, filler.Re, filler.fe);
		filler.insert(2 * enodes->size());

		iterator.temperature        += enodes->size();
		iterator.initialTemperature += enodes->size();
		iterator.coordinates        += enodes->size() * 2;
		iterator.acceleration       += enodes->size() * 2;
		iterator.angularVelocity    += enodes->size() * 3;
		iterator.thickness          += enodes->size();
	}
}

void StructuralMechanics2DController::processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler)
{
	if (info::mesh->boundaryRegions[rindex]->dimension != 1) {
		return;
	}

	auto enodes = info::mesh->boundaryRegions[rindex]->procNodes->cbegin() + filler.begin;
	StructuralMechanics2DKernel::BoundaryIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->boundaryRegions[rindex]->procNodes->datatarray().begin();
	iterator.coordinates = _boundaries[rindex].coordinate.data->datatarray().begin() + noffset * 2;
	iterator.thickness   = _boundaries[rindex].thickness.data->datatarray().begin() + noffset;

	iterator.normalPressure = _boundaries[rindex].normalPressure.data ? _boundaries[rindex].normalPressure.data->datatarray().begin() + noffset : NULL;

	for (esint e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = info::mesh->boundaryRegions[rindex]->epointers->datatarray()[e];

		_kernel->processEdge(matrices, parameters, iterator, filler.Ke, filler.fe);
		filler.insert(2 * enodes->size());

		iterator.coordinates += enodes->size() * 2;
		iterator.thickness   += enodes->size();
		if (iterator.normalPressure) {
			iterator.normalPressure += enodes->size();
		}
	}
}

void StructuralMechanics2DController::processSolution()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		auto enodes = info::mesh->elements->procNodes->cbegin(t);
		StructuralMechanics2DKernel::SolutionIterator iterator;

		iterator.temperature = _ktemperature.data->datatarray().begin(t);
		iterator.coordinates = _kcoordinate.data->datatarray().begin(t);
		iterator.thickness   = _kthickness.data->datatarray().begin(t);

		for (size_t e = info::mesh->elements->distribution[t]; e < info::mesh->elements->distribution[t + 1]; ++e, ++enodes) {
			iterator.element = info::mesh->elements->epointers->datatarray()[e];
			iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

			_kernel->processSolution(iterator);

			iterator.temperature += enodes->size();
			iterator.coordinates += enodes->size() * 2;
			iterator.thickness   += enodes->size();
		}
	}
}

