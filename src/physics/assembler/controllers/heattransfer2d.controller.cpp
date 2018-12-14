
#include "heattransfer2d.controller.h"
#include "../kernels/heattransfer2d.kernel.h"

#include "../../../globals/run.h"
#include "../../../globals/time.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../config/ecf/root.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/nodestore.h"
#include "../../dataholder.h"

using namespace espreso;

HeatTransfer2DControler::HeatTransfer2DControler(HeatTransferLoadStepConfiguration &configuration)
: HeatTransferControler(configuration)
{
	_kernel = new HeatTransfer2DKernel();

	Point defaultMotion(0, 0, 0);
	double defaultHeat = 0; //273.15;
	double defaultThickness = 1;

	_ntemperature.data = new serializededata<eslocal, double>(1, _nDistribution);
	_ncoordinate.data = new serializededata<eslocal, double>(2, _nDistribution);

	_nmotion.isConts = setDefault(configuration.translation_motions, defaultMotion) && defaultMotion.x == defaultMotion.y;
	_nmotion.data = new serializededata<eslocal, double>(2, _nDistribution, defaultMotion.x);

	_nheat.isConts = setDefault(configuration.heat_source, defaultHeat);
	_nheat.data = new serializededata<eslocal, double>(1, _nDistribution, defaultHeat);

	_nthickness.isConts = setDefault(run::ecf->heat_transfer_2d.thickness, defaultThickness);
	_nthickness.data = new serializededata<eslocal, double>(1, _nDistribution, defaultThickness);

	_temperature = run::mesh->nodes->appendData(1, { "TEMPERATURE" });
	_avgThickness = run::mesh->nodes->appendData(1, { }); // printed on elements
	if (run::mesh->hasPhaseChange()) {
		_phaseChange = run::mesh->nodes->appendData(1, { "PHASE" });
		_latentHeat = run::mesh->nodes->appendData(1, { "LATENT_HEAT" });
	}

	if (run::ecf->output.results_selection.translation_motions && configuration.translation_motions.size()) {
		_motion = run::mesh->elements->appendData(2, { "TRANSLATION_MOTION", "TRANSLATION_MOTION_X", "TRANSLATION_MOTION_Y" });
	}

	if (run::ecf->output.results_selection.gradient || run::ecf->heat_transfer_2d.diffusion_split) {
		_gradient = run::mesh->elements->appendData(2, { "GRADIENT", "GRADIENT_X", "GRADIENT_Y" });
	}

	if (run::ecf->output.results_selection.flux) {
		_flux = run::mesh->elements->appendData(2, { "FLUX", "FLUX_X", "FLUX_Y" });
	}

	_boundaries.resize(run::mesh->boundaryRegions.size());
}

HeatTransfer2DControler::~HeatTransfer2DControler()
{
	delete _kernel;
}

void HeatTransfer2DControler::initData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto c = _ncoordinate.data->begin(t);
		for (auto n = run::mesh->elements->procNodes->datatarray().begin(t); n != run::mesh->elements->procNodes->datatarray().end(t); ++n, ++c) {
			c->at(0) = run::mesh->nodes->coordinates->datatarray()[*n].x;
			c->at(1) = run::mesh->nodes->coordinates->datatarray()[*n].y;
		}
	}

	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = time::current;

	updateERegions(run::ecf->heat_transfer_2d.initial_temperature, _ntemperature.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(run::ecf->heat_transfer_2d.thickness, _nthickness.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(_configuration.heat_source, _nheat.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(_configuration.translation_motions, _nmotion.data->datatarray(), 2, cbegin, tbegin, time);

	averageNodeInitilization(_ntemperature.data->datatarray(), _temperature->data);
	averageNodeInitilization(_nthickness.data->datatarray(), _avgThickness->data);

	if (_motion != NULL) {
		nodeValuesToElements(_nmotion.data->datatarray(), _motion->data);
	}

	for (size_t r = 0; r < run::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = run::mesh->boundaryRegions[r];
		if (region->dimension == 1) {

			auto &distribution = region->procNodes->datatarray().distribution();

			_boundaries[r].coordinate.data = new serializededata<eslocal, double>(2, distribution);
			_boundaries[r].temperature.data = new serializededata<eslocal, double>(1, distribution);
			_boundaries[r].thickness.data = new serializededata<eslocal, double>(1, distribution);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto c = _boundaries[r].coordinate.data->begin(t);
				auto temp = _boundaries[r].temperature.data->begin(t);
				auto thick = _boundaries[r].thickness.data->begin(t);
				for (auto n = region->procNodes->datatarray().begin(t); n != region->procNodes->datatarray().end(t); ++n, ++c, ++temp, ++thick) {
					c->at(0) = run::mesh->nodes->coordinates->datatarray()[*n].x;
					c->at(1) = run::mesh->nodes->coordinates->datatarray()[*n].y;
					temp->at(0) = _temperature->data[*n];
					thick->at(0) = _avgThickness->data[*n];
				}
			}

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();
			tbegin = _boundaries[r].temperature.data->datatarray().begin();

			auto flow = _configuration.heat_flow.find(region->name);
			if (flow != _configuration.heat_flow.end()) {
				updateBRegions(flow->second, _boundaries[r].heatFlow, distribution, 2, cbegin, tbegin, time);
				_boundaries[r].regionArea = run::mesh->boundaryRegions[r]->area;
			}
			auto flux = _configuration.heat_flux.find(region->name);
			if (flux != _configuration.heat_flux.end()) {
				updateBRegions(flow->second, _boundaries[r].heatFlux, distribution, 2, cbegin, tbegin, time);
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
				_boundaries[r].htc.data = new serializededata<eslocal, double>(1, distribution);

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

void HeatTransfer2DControler::nextTime()
{
	if (time::isInitial()) {
		return;
	}

	parametersChanged();
}

void HeatTransfer2DControler::parametersChanged()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto temp = _ntemperature.data->datatarray().begin(t);
		for (auto n = run::mesh->elements->procNodes->datatarray().cbegin(t); n != run::mesh->elements->procNodes->datatarray().cend(t); ++n, ++temp) {
			*temp = _temperature->data[*n];
		}
	}

	double *cbegin = _ncoordinate.data->datatarray().data();
	double *tbegin = NULL;
	double time = time::current;

	updateERegions(run::ecf->heat_transfer_2d.thickness, _nthickness.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(_configuration.heat_source, _nheat.data->datatarray(), 2, cbegin, tbegin, time);
	updateERegions(_configuration.translation_motions, _nmotion.data->datatarray(), 2, cbegin, tbegin, time);

	averageNodeInitilization(_nthickness.data->datatarray(), _avgThickness->data);

	if (_motion != NULL) {
		nodeValuesToElements(_nmotion.data->datatarray(), _motion->data);
	}

	for (size_t r = 0; r < run::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = run::mesh->boundaryRegions[r];
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
				updateBRegions(flow->second, _boundaries[r].heatFlux, distribution, 2, cbegin, tbegin, time);
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
				_boundaries[r].htc.data = new serializededata<eslocal, double>(1, distribution);

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

void HeatTransfer2DControler::processElements(Matrices matrices, InstanceFiller &filler)
{
	auto enodes = run::mesh->elements->procNodes->cbegin() + filler.begin;
	HeatTransfer2DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - run::mesh->elements->procNodes->datatarray().begin();
	iterator.temperature = _ntemperature.data->datatarray().begin() + noffset;
	iterator.coordinates = _ncoordinate.data->datatarray().begin() + noffset * 2;
	iterator.motion      = _nmotion.data->datatarray().begin() + noffset * 2;
	iterator.heat        = _nheat.data->datatarray().begin() + noffset;
	iterator.thickness   = _nthickness.data->datatarray().begin() + noffset;

	for (eslocal e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = run::mesh->elements->epointers->datatarray()[e];
		iterator.material = run::mesh->materials[run::mesh->elements->material->datatarray()[e]];

		_kernel->processElement(matrices, iterator, filler.Ke, filler.Me, filler.Re, filler.fe);
		filler.insert(enodes->size());

		iterator.temperature += enodes->size();
		iterator.coordinates += enodes->size() * 2;
		iterator.motion      += enodes->size() * 2;
		iterator.heat        += enodes->size();
		iterator.thickness   += enodes->size();
	}
}

void HeatTransfer2DControler::processBoundary(Matrices matrices, size_t rindex, InstanceFiller &filler)
{
	if (run::mesh->boundaryRegions[rindex]->dimension != 1) {
		return;
	}

	auto enodes = run::mesh->boundaryRegions[rindex]->procNodes->cbegin() + filler.begin;
	HeatTransfer2DKernel::BoundaryIterator iterator;

	size_t noffset = enodes->begin() - run::mesh->boundaryRegions[rindex]->procNodes->datatarray().begin();
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

	for (eslocal e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element = run::mesh->boundaryRegions[rindex]->epointers->datatarray()[e];

		_kernel->processEdge(matrices, iterator, filler.Ke, filler.fe);
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

void HeatTransfer2DControler::processSolution()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		auto enodes = run::mesh->elements->procNodes->cbegin(t);
		HeatTransfer2DKernel::SolutionIterator iterator;

		size_t noffset = enodes->begin() - run::mesh->elements->procNodes->datatarray().begin(t);
		iterator.temperature = _ntemperature.data->datatarray().begin(t);
		iterator.coordinates = _ncoordinate.data->datatarray().begin(t);
		iterator.motion      = _nmotion.data->datatarray().begin(t);
		iterator.heat        = _nheat.data->datatarray().begin(t);
		iterator.thickness   = _nthickness.data->datatarray().begin(t);

		if (run::mesh->hasPhaseChange()) {
			iterator.phase = _phaseChange->data.data() + noffset;
			iterator.latentHeat = _latentHeat->data.data() + noffset;
		}

		if (run::ecf->output.results_selection.gradient) {
			iterator.gradient = _gradient->data.data() + run::mesh->elements->distribution[t] * 2;
		}
		if (run::ecf->output.results_selection.flux) {
			iterator.flux = _flux->data.data() + run::mesh->elements->distribution[t] * 2;
		}

		for (size_t e = run::mesh->elements->distribution[t]; e < run::mesh->elements->distribution[t + 1]; ++e, ++enodes) {
			iterator.element = run::mesh->elements->epointers->datatarray()[e];
			iterator.material = run::mesh->materials[run::mesh->elements->material->datatarray()[e]];

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

			if (run::ecf->output.results_selection.gradient) {
				iterator.gradient += 2;
			}
			if (run::ecf->output.results_selection.flux) {
				iterator.flux += 2;
			}
		}
	}
}

