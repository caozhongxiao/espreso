
#include "esinfo/timeinfo.h"
#include "physics/assembler/dataholder.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "heattransfer2d.controller.h"
#include "physics/assembler/kernels/heattransfer2d.kernel.h"

#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "config/ecf/root.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

HeatTransfer2DController::HeatTransfer2DController(HeatTransfer2DController* previous, HeatTransferLoadStepConfiguration &configuration)
: HeatTransferController(2, previous, configuration), _kernel(new HeatTransfer2DKernel())
{
	if (previous) {
		takeKernelParam(_kcoordinate, previous->_kcoordinate);
		takeKernelParam(_ktemperature, previous->_ktemperature);
		takeKernelParam(_kthickness, previous->_kthickness);
	} else {
		_ntemperature = info::mesh->nodes->appendData(1, { "TEMPERATURE" });
		if (info::mesh->hasPhaseChange()) {
			_ephaseChange = info::mesh->elements->appendData(1, { "PHASE" });
			_elatentHeat = info::mesh->elements->appendData(1, { "LATENT_HEAT" });
		}

		if (info::ecf->output.results_selection.thickness) {
			_ethickness = info::mesh->elements->appendData(1, { "THICKNESS" });
		}

		if (info::ecf->output.results_selection.translation_motions) {
			for (auto it = info::ecf->heat_transfer_2d.load_steps_settings.begin(); it != info::ecf->heat_transfer_2d.load_steps_settings.end(); ++it) {
				if (it->second.translation_motions.size()) {
					_emotion = info::mesh->elements->appendData(2, { "TRANSLATION_MOTION", "TRANSLATION_MOTION_X", "TRANSLATION_MOTION_Y" });
					break;
				}
			}
		}

		if (info::ecf->output.results_selection.gradient || info::ecf->heat_transfer_2d.diffusion_split) {
			_egradient = info::mesh->elements->appendData(2, { "GRADIENT", "GRADIENT_X", "GRADIENT_Y" });
		}

		if (info::ecf->output.results_selection.flux) {
			_eflux = info::mesh->elements->appendData(2, { "FLUX", "FLUX_X", "FLUX_Y" });
		}

		setCoordinates(_kcoordinate);
		initKernelParam(_ktemperature, info::ecf->heat_transfer_2d.initial_temperature, 0);
		initKernelParam(_kthickness, info::ecf->heat_transfer_2d.thickness, 1);

		// initial temperature is defined for element regions -> set for elements, average and set dirichlet
		if (!_ktemperature.isConts) {
			updateKernelParam(_ktemperature, info::ecf->heat_transfer_2d.initial_temperature, _kcoordinate.data->datatarray().data(), NULL);
		}
		kernelToNodes(_ktemperature, _ntemperature->data);
		if (info::ecf->heat_transfer_2d.init_temp_respect_bc) {
			setDirichlet(_ntemperature->data);
		}

		if (_ethickness && _kthickness.isConts) {
			kernelToElements(_kthickness, _ethickness->data);
		}
	}

	initKernelParam(_kmotion, configuration.translation_motions, 0);
	initKernelParam(_kheat, configuration.heat_source, 0);

	if (_emotion && _kmotion.isConts) {
		kernelToElements(_kmotion, _emotion->data);
	}

	_boundaries.resize(info::mesh->boundaryRegions.size(), BoundaryParameters(2));
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 1) {
			_boundaries[r].regionArea = info::mesh->boundaryRegions[r]->area;
			setCoordinates(_boundaries[r].coordinate, region->procNodes);
			initKernelParam(_boundaries[r].temperature, info::ecf->heat_transfer_2d.initial_temperature, 0, region);
			initKernelParam(_boundaries[r].thickness, info::ecf->heat_transfer_2d.thickness, 1, region);

			auto heatFlow = _configuration.heat_flow.find(region->name);
			if (heatFlow != _configuration.heat_flow.end()) {
				initKernelParam(_boundaries[r].heatFlow, heatFlow->second, 0, region);
			}
			auto heatFlux = _configuration.heat_flux.find(region->name);
			if (heatFlux != _configuration.heat_flux.end()) {
				initKernelParam(_boundaries[r].heatFlux, heatFlux->second, 0, region);
			}

			auto radiation = _configuration.diffuse_radiation.find(region->name);
			if (radiation != _configuration.diffuse_radiation.end()) {
				initKernelParam(_boundaries[r].emissivity, radiation->second.emissivity, 0, region);
				initKernelParam(_boundaries[r].externalTemperature, radiation->second.external_temperature, 0, region);
			}

			auto convection = _configuration.convection.find(region->name);
			if (convection != _configuration.convection.end()) {
				if (_boundaries[r].externalTemperature.data == NULL) {
					initKernelParam(_boundaries[r].externalTemperature, convection->second.external_temperature, 0, region);
				}
				initKernelParam(_boundaries[r].htc, convection->second.heat_transfer_coefficient, 0, region);
			}
		}
	}
}

HeatTransfer2DController::~HeatTransfer2DController()
{
	if (_kernel) {
		delete _kernel;
	}
}

const PhysicsConfiguration& HeatTransfer2DController::configuration() const
{
	return info::ecf->heat_transfer_2d;
}

void HeatTransfer2DController::nextTime()
{
	parametersChanged();
}

void HeatTransfer2DController::parametersChanged()
{
	nodesToKernels(_ntemperature->data, _ktemperature);

	double *cbegin = _kcoordinate.data->datatarray().data();
	double *tbegin = _ktemperature.data->datatarray().data();

	if (!_kthickness.isConts) {
		updateKernelParam(_kthickness, info::ecf->heat_transfer_2d.thickness, cbegin, tbegin);
		if (_ethickness) {
			kernelToElements(_kthickness, _ethickness->data);
		}
	}
	if (!_kheat.isConts) {
		updateKernelParam(_kheat, _configuration.heat_source, cbegin, tbegin);
	}
	if (!_kmotion.isConts) {
		updateKernelParam(_kmotion, _configuration.translation_motions, cbegin, tbegin);
		if (_emotion) {
			kernelToElements(_kmotion, _emotion->data);
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 1) {

			nodesToKernels(_ntemperature->data, _boundaries[r].temperature, region->procNodes);

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();
			tbegin = _boundaries[r].temperature.data->datatarray().begin();

			if (!_boundaries[r].thickness.isConts) {
				kernelToBoundary(_kthickness, _boundaries[r].thickness, region);
			}

			auto heatFlow = _configuration.heat_flow.find(region->name);
			if (heatFlow != _configuration.heat_flow.end()) {
				if (!_boundaries[r].heatFlow.isConts) {
					updateKernelParam(_boundaries[r].heatFlow, heatFlow->second, cbegin, tbegin, region);
				}
			}
			auto heatFlux = _configuration.heat_flux.find(region->name);
			if (heatFlux != _configuration.heat_flux.end()) {
				if (!_boundaries[r].heatFlux.isConts) {
					updateKernelParam(_boundaries[r].heatFlux, heatFlux->second, cbegin, tbegin, region);
				}
			}

			auto radiation = _configuration.diffuse_radiation.find(region->name);
			if (radiation != _configuration.diffuse_radiation.end()) {
				if (!_boundaries[r].emissivity.isConts) {
					updateKernelParam(_boundaries[r].emissivity, radiation->second.emissivity, cbegin, tbegin, region);
				}
				if (!_boundaries[r].externalTemperature.isConts) {
					updateKernelParam(_boundaries[r].externalTemperature, radiation->second.external_temperature, cbegin, tbegin, region);
				}
			}

			auto convection = _configuration.convection.find(region->name);
			if (convection != _configuration.convection.end()) {
				if (!_boundaries[r].externalTemperature.isConts) {
					updateKernelParam(_boundaries[r].externalTemperature, convection->second.external_temperature, cbegin, tbegin, region);
				}
				#pragma omp parallel for
				for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
					for (size_t i = region->distribution[t]; i < region->distribution[t + 1]; ++i) {
						_boundaries[r].htc.data->datatarray()[i] = _kernel->convectionHTC(convection->second, 2, cbegin + i * 2, time::current, *(tbegin + i));
					}
				}
			}
		}
	}

	if (info::ecf->heat_transfer_2d.diffusion_split) {
		processSolution(); // compute gradient
	}
}

void HeatTransfer2DController::processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler)
{
	auto enodes = info::mesh->elements->procNodes->cbegin() + filler.begin;
	HeatTransfer2DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin();
	iterator.temperature = _ktemperature.data->datatarray().begin() + noffset;
	iterator.coordinates = _kcoordinate.data->datatarray().begin() + noffset * 2;
	iterator.motion      = _kmotion.data->datatarray().begin() + noffset * 2;
	iterator.heat        = _kheat.data->datatarray().begin() + noffset;
	iterator.thickness   = _kthickness.data->datatarray().begin() + noffset;

	for (esint e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = info::mesh->elements->epointers->datatarray()[e];
		iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

		_kernel->processElement(matrices, parameters, iterator, filler.Ke, filler.Me, filler.Re, filler.fe);
		filler.insert(enodes->size());

		iterator.temperature += enodes->size();
		iterator.coordinates += enodes->size() * 2;
		iterator.motion      += enodes->size() * 2;
		iterator.heat        += enodes->size();
		iterator.thickness   += enodes->size();
	}
}

void HeatTransfer2DController::processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler)
{
	if (info::mesh->boundaryRegions[rindex]->dimension != 1) {
		return;
	}

	auto enodes = info::mesh->boundaryRegions[rindex]->procNodes->cbegin() + filler.begin;
	HeatTransfer2DKernel::BoundaryIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->boundaryRegions[rindex]->procNodes->datatarray().begin();
	iterator.temperature = _boundaries[rindex].temperature.data->datatarray().begin() + noffset;
	iterator.coordinates = _boundaries[rindex].coordinate.data->datatarray().begin() + noffset * 2;
	iterator.thickness   = _boundaries[rindex].thickness.data->datatarray().begin() + noffset;

	iterator.regionArea = _boundaries[rindex].regionArea;
	iterator.heatFlow = _boundaries[rindex].heatFlow.data ? _boundaries[rindex].heatFlow.data->datatarray().begin() + noffset : NULL;
	iterator.heatFlux = _boundaries[rindex].heatFlux.data ? _boundaries[rindex].heatFlux.data->datatarray().begin() + noffset : NULL;
	iterator.emissivity = _boundaries[rindex].emissivity.data ? _boundaries[rindex].emissivity.data->datatarray().begin() + noffset : NULL;
	iterator.externalTemperature = _boundaries[rindex].externalTemperature.data ? _boundaries[rindex].externalTemperature.data->datatarray().begin() + noffset : NULL;
	iterator.htc = _boundaries[rindex].htc.data ? _boundaries[rindex].htc.data->datatarray().begin() + noffset : NULL;
	iterator.radiation = iterator.emissivity != NULL;
	iterator.convection = iterator.htc != NULL;

	for (esint e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = info::mesh->boundaryRegions[rindex]->epointers->datatarray()[e];

		_kernel->processEdge(matrices, parameters, iterator, filler.Ke, filler.fe);
		filler.insert(enodes->size());

		iterator.temperature += enodes->size();
		iterator.coordinates += enodes->size() * 2;
		iterator.thickness   += enodes->size();
		if (iterator.heatFlow) {
			iterator.heatFlow += enodes->size();
		}
		if (iterator.heatFlux) {
			iterator.heatFlux += enodes->size();
		}
		if (iterator.emissivity) {
			iterator.emissivity += enodes->size();
		}
		if (iterator.htc) {
			iterator.htc += enodes->size();
		}
		if (iterator.externalTemperature) {
			iterator.externalTemperature += enodes->size();
		}
	}
}

void HeatTransfer2DController::processSolution()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		auto enodes = info::mesh->elements->procNodes->cbegin(t);
		HeatTransfer2DKernel::SolutionIterator iterator;

		iterator.temperature = _ktemperature.data->datatarray().begin(t);
		iterator.coordinates = _kcoordinate.data->datatarray().begin(t);
		iterator.motion      = _kmotion.data->datatarray().begin(t);
		iterator.heat        = _kheat.data->datatarray().begin(t);
		iterator.thickness   = _kthickness.data->datatarray().begin(t);

		if (info::mesh->hasPhaseChange()) {
			iterator.phase = _ephaseChange->data.data() + info::mesh->elements->distribution[t];
			iterator.latentHeat = _elatentHeat->data.data() + info::mesh->elements->distribution[t];
		}

		if (info::ecf->output.results_selection.gradient) {
			iterator.gradient = _egradient->data.data() + info::mesh->elements->distribution[t] * 2;
		}
		if (info::ecf->output.results_selection.flux) {
			iterator.flux = _eflux->data.data() + info::mesh->elements->distribution[t] * 2;
		}

		for (size_t e = info::mesh->elements->distribution[t]; e < info::mesh->elements->distribution[t + 1]; ++e, ++enodes) {
			iterator.element = info::mesh->elements->epointers->datatarray()[e];
			iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

			_kernel->processSolution(iterator);

			iterator.temperature += enodes->size();
			iterator.coordinates += enodes->size() * 2;
			iterator.motion      += enodes->size() * 2;
			iterator.heat        += enodes->size();
			iterator.thickness   += enodes->size();

			if (iterator.phase) {
				iterator.phase += 1;
				iterator.latentHeat += 1;
			}

			if (info::ecf->output.results_selection.gradient) {
				iterator.gradient += 2;
			}
			if (info::ecf->output.results_selection.flux) {
				iterator.flux += 2;
			}
		}
	}
}

