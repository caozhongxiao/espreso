
#include "physics/assembler/dataholder.h"
#include "esinfo/time.h"
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

HeatTransfer2DController::HeatTransfer2DController(HeatTransferLoadStepConfiguration &configuration)
: HeatTransferController(configuration)
{
	_kernel = new HeatTransfer2DKernel();

	Point defaultMotion(0, 0, 0);
	double defaultHeat = 0; //273.15;
	double defaultThickness = 1;

	_ntemperature.data = new serializededata<esint, double>(1, _nDistribution);
	_ncoordinate.data = new serializededata<esint, double>(2, _nDistribution);

	_nmotion.isConts = setDefault(configuration.translation_motions, defaultMotion) && defaultMotion.x == defaultMotion.y;
	_nmotion.data = new serializededata<esint, double>(2, _nDistribution, defaultMotion.x);

	_nheat.isConts = setDefault(configuration.heat_source, defaultHeat);
	_nheat.data = new serializededata<esint, double>(1, _nDistribution, defaultHeat);

	_nthickness.isConts = setDefault(info::ecf->heat_transfer_2d.thickness, defaultThickness);
	_nthickness.data = new serializededata<esint, double>(1, _nDistribution, defaultThickness);

	_temperature = info::mesh->nodes->appendData(1, { "TEMPERATURE" });
	_avgThickness = info::mesh->nodes->appendData(1, { }); // printed on elements
	if (info::mesh->hasPhaseChange()) {
		_phaseChange = info::mesh->elements->appendData(1, { "PHASE" });
		_latentHeat = info::mesh->elements->appendData(1, { "LATENT_HEAT" });
	}

	if (info::ecf->output.results_selection.translation_motions && configuration.translation_motions.size()) {
		_motion = info::mesh->elements->appendData(2, { "TRANSLATION_MOTION", "TRANSLATION_MOTION_X", "TRANSLATION_MOTION_Y" });
	}

	if (info::ecf->output.results_selection.gradient || info::ecf->heat_transfer_2d.diffusion_split) {
		_gradient = info::mesh->elements->appendData(2, { "GRADIENT", "GRADIENT_X", "GRADIENT_Y" });
	}

	if (info::ecf->output.results_selection.flux) {
		_flux = info::mesh->elements->appendData(2, { "FLUX", "FLUX_X", "FLUX_Y" });
	}

	_boundaries.resize(info::mesh->boundaryRegions.size());
}

HeatTransfer2DController::~HeatTransfer2DController()
{
	delete _kernel;
}

const PhysicsConfiguration& HeatTransfer2DController::configuration() const
{
	return info::ecf->heat_transfer_2d;
}

void HeatTransfer2DController::initData()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto c = _ncoordinate.data->begin(t);
		for (auto n = info::mesh->elements->procNodes->datatarray().begin(t); n != info::mesh->elements->procNodes->datatarray().end(t); ++n, ++c) {
			c->at(0) = info::mesh->nodes->coordinates->datatarray()[*n].x;
			c->at(1) = info::mesh->nodes->coordinates->datatarray()[*n].y;
		}
	}

	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = time::current;

	updateERegions(info::ecf->heat_transfer_2d.initial_temperature, _ntemperature.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(info::ecf->heat_transfer_2d.thickness, _nthickness.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(_configuration.heat_source, _nheat.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(_configuration.translation_motions, _nmotion.data->datatarray(), 2, cbegin, tbegin, time);

	if (info::ecf->heat_transfer_2d.init_temp_respect_bc) {
		initDirichletData(_ntemperature.data->datatarray());
	}
	averageNodeInitilization(_ntemperature.data->datatarray(), _temperature->data);
	averageNodeInitilization(_nthickness.data->datatarray(), _avgThickness->data);

	if (_motion != NULL) {
		nodeValuesToElements(2, _nmotion.data->datatarray(), _motion->data);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 1) {

			auto &distribution = region->procNodes->datatarray().distribution();

			_boundaries[r].coordinate.data = new serializededata<esint, double>(2, distribution);
			_boundaries[r].temperature.data = new serializededata<esint, double>(1, distribution);
			_boundaries[r].thickness.data = new serializededata<esint, double>(1, distribution);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto c = _boundaries[r].coordinate.data->begin(t);
				auto temp = _boundaries[r].temperature.data->begin(t);
				auto thick = _boundaries[r].thickness.data->begin(t);
				for (auto n = region->procNodes->datatarray().begin(t); n != region->procNodes->datatarray().end(t); ++n, ++c, ++temp, ++thick) {
					c->at(0) = info::mesh->nodes->coordinates->datatarray()[*n].x;
					c->at(1) = info::mesh->nodes->coordinates->datatarray()[*n].y;
					temp->at(0) = _temperature->data[*n];
					thick->at(0) = _avgThickness->data[*n];
				}
			}

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();
			tbegin = _boundaries[r].temperature.data->datatarray().begin();

			auto flow = _configuration.heat_flow.find(region->name);
			if (flow != _configuration.heat_flow.end()) {
				_boundaries[r].heatFlow.data = new serializededata<esint, double>(1, distribution);
				_boundaries[r].heatFlow.isConts = flow->second.evaluator->isConstant();
				updateBRegions(flow->second, _boundaries[r].heatFlow, distribution, 2, cbegin, tbegin, time);
				_boundaries[r].regionArea = info::mesh->boundaryRegions[r]->area;
			}
			auto flux = _configuration.heat_flux.find(region->name);
			if (flux != _configuration.heat_flux.end()) {
				_boundaries[r].heatFlux.data = new serializededata<esint, double>(1, distribution);
				_boundaries[r].heatFlux.isConts = flux->second.evaluator->isConstant();
				updateBRegions(flux->second, _boundaries[r].heatFlux, distribution, 2, cbegin, tbegin, time);
			}

			auto radiation = _configuration.diffuse_radiation.find(region->name);
			if (radiation != _configuration.diffuse_radiation.end()) {
				_boundaries[r].emissivity.data = new serializededata<esint, double>(1, distribution);
				_boundaries[r].emissivity.isConts = radiation->second.emissivity.evaluator->isConstant();
				_boundaries[r].externalTemperature.data = new serializededata<esint, double>(1, distribution);
				_boundaries[r].externalTemperature.isConts = radiation->second.external_temperature.evaluator->isConstant();
				updateBRegions(radiation->second.emissivity, _boundaries[r].emissivity, distribution, 2, cbegin, tbegin, time);
				updateBRegions(radiation->second.external_temperature, _boundaries[r].externalTemperature, distribution, 2, cbegin, tbegin, time);
			}

			auto convection = _configuration.convection.find(region->name);
			if (convection != _configuration.convection.end()) {
				if (_boundaries[r].externalTemperature.data == NULL) {
					_boundaries[r].externalTemperature.data = new serializededata<esint, double>(1, distribution);
					_boundaries[r].externalTemperature.isConts = convection->second.external_temperature.evaluator->isConstant();
					updateBRegions(convection->second.external_temperature, _boundaries[r].externalTemperature, distribution, 2, cbegin, tbegin, time);
				}
				_boundaries[r].htc.data = new serializededata<esint, double>(1, distribution);

				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
						_boundaries[r].htc.data->datatarray()[i] = _kernel->convectionHTC(convection->second, 2, cbegin + i * 2, time, *(tbegin + i));
					}
				}
			}
		}
	}
}

void HeatTransfer2DController::nextTime()
{
	if (time::isInitial()) {
		return;
	}

	parametersChanged();
}

void HeatTransfer2DController::parametersChanged()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto temp = _ntemperature.data->datatarray().begin(t);
		for (auto n = info::mesh->elements->procNodes->datatarray().cbegin(t); n != info::mesh->elements->procNodes->datatarray().cend(t); ++n, ++temp) {
			*temp = _temperature->data[*n];
		}
	}

	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = time::current;

	updateERegions(info::ecf->heat_transfer_2d.thickness, _nthickness.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(_configuration.heat_source, _nheat.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(_configuration.translation_motions, _nmotion.data->datatarray(), 2, cbegin, tbegin, time);

	averageNodeInitilization(_nthickness.data->datatarray(), _avgThickness->data);

	if (_motion != NULL) {
		nodeValuesToElements(2, _nmotion.data->datatarray(), _motion->data);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 1) {

			auto &distribution = region->procNodes->datatarray().distribution();

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto temp = _boundaries[r].temperature.data->begin(t);
				auto thick = _boundaries[r].thickness.data->begin(t);
				for (auto n = region->procNodes->datatarray().begin(t); n != region->procNodes->datatarray().end(t); ++n, ++temp, ++thick) {
					temp->at(0) = _temperature->data[*n];
					thick->at(0) = _avgThickness->data[*n];
				}
			}

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();
			tbegin = _boundaries[r].temperature.data->datatarray().begin();

			auto flow = _configuration.heat_flow.find(region->name);
			if (flow != _configuration.heat_flow.end()) {
				updateBRegions(flow->second, _boundaries[r].heatFlow, distribution, 2, cbegin, tbegin, time);
			}
			auto flux = _configuration.heat_flux.find(region->name);
			if (flux != _configuration.heat_flux.end()) {
				updateBRegions(flux->second, _boundaries[r].heatFlux, distribution, 2, cbegin, tbegin, time);
			}

			auto radiation = _configuration.diffuse_radiation.find(region->name);
			if (radiation != _configuration.diffuse_radiation.end()) {
				updateBRegions(radiation->second.emissivity, _boundaries[r].emissivity, distribution, 2, cbegin, tbegin, time);
				updateBRegions(radiation->second.external_temperature, _boundaries[r].externalTemperature, distribution, 2, cbegin, tbegin, time);
			}

			auto convection = _configuration.convection.find(region->name);
			if (convection != _configuration.convection.end()) {
				if (_boundaries[r].externalTemperature.data == NULL) {
					updateBRegions(convection->second.external_temperature, _boundaries[r].externalTemperature, distribution, 2, cbegin, tbegin, time);
				}
				_boundaries[r].htc.data = new serializededata<esint, double>(1, distribution);

				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
						_boundaries[r].htc.data->datatarray()[i] = _kernel->convectionHTC(convection->second, 2, cbegin + i * 2, time, *(tbegin + i));
					}
				}
			}
		}
	}
}

void HeatTransfer2DController::processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler)
{
	auto enodes = info::mesh->elements->procNodes->cbegin() + filler.begin;
	HeatTransfer2DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin();
	iterator.temperature = _ntemperature.data->datatarray().begin() + noffset;
	iterator.coordinates = _ncoordinate.data->datatarray().begin() + noffset * 2;
	iterator.motion      = _nmotion.data->datatarray().begin() + noffset * 2;
	iterator.heat        = _nheat.data->datatarray().begin() + noffset;
	iterator.thickness   = _nthickness.data->datatarray().begin() + noffset;

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

		iterator.temperature = _ntemperature.data->datatarray().begin(t);
		iterator.coordinates = _ncoordinate.data->datatarray().begin(t);
		iterator.motion      = _nmotion.data->datatarray().begin(t);
		iterator.heat        = _nheat.data->datatarray().begin(t);
		iterator.thickness   = _nthickness.data->datatarray().begin(t);

		if (info::mesh->hasPhaseChange()) {
			iterator.phase = _phaseChange->data.data() + info::mesh->elements->distribution[t];
			iterator.latentHeat = _latentHeat->data.data() + info::mesh->elements->distribution[t];
		}

		if (info::ecf->output.results_selection.gradient) {
			iterator.gradient = _gradient->data.data() + info::mesh->elements->distribution[t] * 2;
		}
		if (info::ecf->output.results_selection.flux) {
			iterator.flux = _flux->data.data() + info::mesh->elements->distribution[t] * 2;
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

