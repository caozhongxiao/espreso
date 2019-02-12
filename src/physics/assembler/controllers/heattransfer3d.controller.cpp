
#include "heattransfer3d.controller.h"

#include "esinfo/time.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/kernels/heattransfer3d.kernel.h"

#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/surfacestore.h"

#include "wrappers/bem/bemwrapper.h"

using namespace espreso;

HeatTransfer3DController::HeatTransfer3DController(HeatTransfer3DController* previous, HeatTransferLoadStepConfiguration &configuration)
: HeatTransferController(3, previous, configuration), _kernel(new HeatTransfer3DKernel()), _bem(info::mesh->elements->ndomains, NULL)
{
	// TODO: optimize controller when only BEM is needed

	if (previous) {
		takeKernelParam(_kcoordinate, previous->_kcoordinate);
		takeKernelParam(_ktemperature, previous->_ktemperature);
	} else {
		_ntemperature = info::mesh->nodes->appendData(1, { "TEMPERATURE" });
		if (info::mesh->hasPhaseChange()) {
			_ephaseChange = info::mesh->elements->appendData(1, { "PHASE" });
			_elatentHeat = info::mesh->elements->appendData(1, { "LATENT_HEAT" });
		}

		if (info::ecf->output.results_selection.translation_motions) {
			for (auto it = info::ecf->heat_transfer_3d.load_steps_settings.begin(); it != info::ecf->heat_transfer_3d.load_steps_settings.end(); ++it) {
				if (it->second.translation_motions.size()) {
					_emotion = info::mesh->elements->appendData(3, { "TRANSLATION_MOTION", "TRANSLATION_MOTION_X", "TRANSLATION_MOTION_Y", "TRANSLATION_MOTION_Z" });
					break;
				}
			}
		}

		if (info::ecf->output.results_selection.gradient || info::ecf->heat_transfer_3d.diffusion_split) {
			_egradient = info::mesh->elements->appendData(3, { "GRADIENT", "GRADIENT_X", "GRADIENT_Y", "GRADIENT_Z" });
		}

		if (info::ecf->output.results_selection.flux) {
			_eflux = info::mesh->elements->appendData(3, { "FLUX", "FLUX_X", "FLUX_Y", "FLUX_Z" });
		}

		setCoordinates(_kcoordinate);
		initKernelParam(_ktemperature, info::ecf->heat_transfer_3d.initial_temperature, 0);

		// initial temperature is defined for element regions -> set for elements, average and set dirichlet
		if (!_ktemperature.isConts) {
			updateKernelParam(_ktemperature, info::ecf->heat_transfer_3d.initial_temperature, _kcoordinate.data->datatarray().data(), NULL);
		}
		kernelToNodes(_ktemperature, _ntemperature->data);
		if (info::ecf->heat_transfer_3d.init_temp_respect_bc) {
			setDirichlet(_ntemperature->data);
		}
	}

	initKernelParam(_kmotion, configuration.translation_motions, 0);
	initKernelParam(_kheat, configuration.heat_source, 0);

	_boundaries.resize(info::mesh->boundaryRegions.size(), BoundaryParameters(3));
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 2) {
			_boundaries[r].regionArea = info::mesh->boundaryRegions[r]->area;
			setCoordinates(_boundaries[r].coordinate, region->procNodes);
			initKernelParam(_boundaries[r].temperature, info::ecf->heat_transfer_3d.initial_temperature, 0, region);

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

HeatTransfer3DController::~HeatTransfer3DController()
{
	delete _kernel;
	for (size_t i = 0; i < _bem.size(); i++) {
		if (_bem[i] != NULL) {
			BEM4I::deleteData(_bem[i]);
		}
	}
}

const PhysicsConfiguration& HeatTransfer3DController::configuration() const
{
	return info::ecf->heat_transfer_3d;
}

void HeatTransfer3DController::nextTime()
{
	parametersChanged();
}

void HeatTransfer3DController::parametersChanged()
{
	nodesToKernels(_ntemperature->data, _ktemperature);

	double *cbegin = _kcoordinate.data->datatarray().data();
	double *tbegin = _ktemperature.data->datatarray().data();

	if (!_kheat.isConts) {
		updateKernelParam(_kheat, _configuration.heat_source, cbegin, tbegin);
	}
	if (!_kmotion.isConts) {
		updateKernelParam(_kmotion, _configuration.translation_motions, cbegin, tbegin);
		kernelToElements(_kmotion, _emotion->data);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 2) {

			nodesToKernels(_ntemperature->data, _boundaries[r].temperature, region->procNodes);

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();
			tbegin = _boundaries[r].temperature.data->datatarray().begin();

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
						_boundaries[r].htc.data->datatarray()[i] = _kernel->convectionHTC(convection->second, 3, cbegin + i * 3, time::current, *(tbegin + i));
					}
				}
			}
		}
	}

	if (info::ecf->heat_transfer_3d.diffusion_split) {
		processSolution(); // compute gradient
	}
}

void HeatTransfer3DController::processBEMdomain(esint domain, double *values)
{
	esint nodes = info::mesh->domainsSurface->cdistribution[domain + 1] - info::mesh->domainsSurface->cdistribution[domain];
	esint noffset = info::mesh->domainsSurface->cdistribution[domain];
	esint elements = info::mesh->domainsSurface->tdistribution[domain + 1] - info::mesh->domainsSurface->tdistribution[domain];
	esint eoffset = info::mesh->domainsSurface->tdistribution[domain];
	auto mat = info::mesh->materials[info::mesh->elements->material->datatarray()[info::mesh->elements->elementsDistribution[domain]]];
	std::vector<double> K(nodes * nodes);

	BEM4I::getLaplace(
			_bem[domain], K.data(),
			nodes, reinterpret_cast<double*>((info::mesh->domainsSurface->coordinates->begin() + noffset)->begin()),
			elements, (info::mesh->domainsSurface->triangles->begin() + eoffset)->begin(),
			mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(Point(0, 0, 0), time::current, 0));

	for (esint r = 0, i = 0; r < nodes; r++) {
		for (esint c = r; c < nodes; c++, i++) {
			values[i] = K[r * nodes + c];
		}
	}
}

void HeatTransfer3DController::fillBEMInterior(esint domain, double *values)
{
	esint nodes = info::mesh->domainsSurface->cdistribution[domain + 1] - info::mesh->domainsSurface->cdistribution[domain];
	esint noffset = info::mesh->domainsSurface->cdistribution[domain];
	esint points = info::mesh->nodes->dintervals[domain].back().end - info::mesh->nodes->dintervals[domain].back().begin;
	esint poffset = info::mesh->nodes->dintervals[domain].back().begin;
	esint elements = info::mesh->domainsSurface->tdistribution[domain + 1] - info::mesh->domainsSurface->tdistribution[domain];
	esint eoffset = info::mesh->domainsSurface->tdistribution[domain];
	auto mat = info::mesh->materials[info::mesh->elements->material->datatarray()[info::mesh->elements->elementsDistribution[domain]]];

	if (nodes == info::mesh->nodes->dintervals[domain].back().DOFOffset + points) {
		return; // no interior nodes
	}

	BEM4I::evaluateLaplace(
			_bem[domain], values + nodes,
			nodes, reinterpret_cast<double*>((info::mesh->domainsSurface->coordinates->begin() + noffset)->begin()),
			elements, (info::mesh->domainsSurface->triangles->begin() + eoffset)->begin(),
			points, reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data() + poffset),
			mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(Point(0, 0, 0), time::current, 0),
			values);
}

void HeatTransfer3DController::processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler)
{
	auto enodes = info::mesh->elements->procNodes->cbegin() + filler.begin;
	HeatTransfer3DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin();
	iterator.temperature = _ktemperature.data->datatarray().begin() + noffset;
	iterator.coordinates = _kcoordinate.data->datatarray().begin() + noffset * 3;
	iterator.motion      = _kmotion.data->datatarray().begin() + noffset * 3;
	iterator.heat        = _kheat.data->datatarray().begin() + noffset;

	for (esint e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = info::mesh->elements->epointers->datatarray()[e];
		iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

		_kernel->processElement(matrices, parameters, iterator, filler.Ke, filler.Me, filler.Re, filler.fe);
		filler.insert(enodes->size());

		iterator.temperature += enodes->size();
		iterator.coordinates += enodes->size() * 3;
		iterator.motion      += enodes->size() * 3;
		iterator.heat        += enodes->size();
	}
}

void HeatTransfer3DController::processBoundary(Matrices matrices, const SolverParameters &parameters, size_t rindex, InstanceFiller &filler)
{
	if (info::mesh->boundaryRegions[rindex]->dimension != 2) {
		return;
	}

	auto enodes = info::mesh->boundaryRegions[rindex]->procNodes->cbegin() + filler.begin;
	HeatTransfer3DKernel::BoundaryIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->boundaryRegions[rindex]->procNodes->datatarray().begin();
	iterator.temperature = _boundaries[rindex].temperature.data->datatarray().begin() + noffset;
	iterator.coordinates = _boundaries[rindex].coordinate.data->datatarray().begin() + noffset * 3;

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

		_kernel->processFace(matrices, parameters, iterator, filler.Ke, filler.fe);
		filler.insert(enodes->size());

		iterator.temperature += enodes->size();
		iterator.coordinates += enodes->size() * 3;
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

void HeatTransfer3DController::processSolution()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		auto enodes = info::mesh->elements->procNodes->cbegin(t);
		HeatTransfer3DKernel::SolutionIterator iterator;

		iterator.temperature = _ktemperature.data->datatarray().begin(t);
		iterator.coordinates = _kcoordinate.data->datatarray().begin(t);
		iterator.motion      = _kmotion.data->datatarray().begin(t);
		iterator.heat        = _kheat.data->datatarray().begin(t);

		if (info::mesh->hasPhaseChange()) {
			iterator.phase = _ephaseChange->data.data() + info::mesh->elements->distribution[t];
			iterator.latentHeat = _elatentHeat->data.data() + info::mesh->elements->distribution[t];
		}

		if (info::ecf->output.results_selection.gradient) {
			iterator.gradient = _egradient->data.data() + info::mesh->elements->distribution[t] * 3;
		}
		if (info::ecf->output.results_selection.flux) {
			iterator.flux = _eflux->data.data() + info::mesh->elements->distribution[t] * 3;
		}

		for (size_t e = info::mesh->elements->distribution[t]; e < info::mesh->elements->distribution[t + 1]; ++e, ++enodes) {
			iterator.element = info::mesh->elements->epointers->datatarray()[e];
			iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

			_kernel->processSolution(iterator);

			iterator.temperature += enodes->size();
			iterator.coordinates += enodes->size() * 3;
			iterator.motion      += enodes->size() * 3;
			iterator.heat        += enodes->size();

			if (iterator.phase) {
				iterator.phase += 1;
				iterator.latentHeat += 1;
			}

			if (info::ecf->output.results_selection.gradient) {
				iterator.gradient += 3;
			}
			if (info::ecf->output.results_selection.flux) {
				iterator.flux += 3;
			}
		}
	}

//	if (BEM) {
//		if (_instance->primalSolution[domain].size() < (*_temperature->decomposedData)[domain].size()) {
//			bem4i::evaluateLaplaceRepresentationFormula(
//					_instance->K[domain].rows,
//					reinterpret_cast<double*>(_mesh->domainsSurface->coordinates->datatarray().data() + _mesh->domainsSurface->cdistribution[domain]),
//					_mesh->domainsSurface->tdistribution[domain + 1] - _mesh->domainsSurface->tdistribution[domain],
//					_mesh->domainsSurface->triangles->datatarray().data() + 3 * _mesh->domainsSurface->tdistribution[domain],
//					_mesh->procNodes->dintervals[domain].back().end - _mesh->procNodes->dintervals[domain].back().begin,
//					reinterpret_cast<double*>(_mesh->procNodes->coordinates->datatarray().data() + _mesh->procNodes->dintervals[domain].back().begin),
//					(*_temperature->decomposedData)[domain].data() + _mesh->procNodes->dintervals[domain].back().DOFOffset,
//					(*_temperature->decomposedData)[domain].data(),
//					4,
//					_BEMData[domain],
//					0);
//		}
//	}
}







