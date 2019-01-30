
#include "physics/assembler/dataholder.h"
#include "esinfo/time.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "heattransfer3d.controller.h"
#include "physics/assembler/kernels/heattransfer3d.kernel.h"

#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

HeatTransfer3DController::HeatTransfer3DController(HeatTransferLoadStepConfiguration &configuration)
: HeatTransferController(configuration)
{
	_kernel = new HeatTransfer3DKernel();

	Point defaultMotion(0, 0, 0);
	double defaultHeat = 0; //273.15;

	_ntemperature.data = new serializededata<esint, double>(1, _nDistribution);
	_ncoordinate.data = new serializededata<esint, double>(3, _nDistribution);

	_nmotion.isConts = setDefault(_configuration.translation_motions, defaultMotion) && defaultMotion.x == defaultMotion.y && defaultMotion.x == defaultMotion.z;
	_nmotion.data = new serializededata<esint, double>(3, _nDistribution, defaultMotion.x);

	_nheat.isConts = setDefault(_configuration.heat_source, defaultHeat);
	_nheat.data = new serializededata<esint, double>(1, _nDistribution, defaultHeat);

	_temperature = info::mesh->nodes->appendData(1, { "TEMPERATURE" });
	if (info::mesh->hasPhaseChange()) {
		_phaseChange = info::mesh->elements->appendData(1, { "PHASE" });
		_latentHeat = info::mesh->elements->appendData(1, { "LATENT_HEAT" });
	}

	if (info::ecf->output.results_selection.translation_motions && _configuration.translation_motions.size()) {
		_motion = info::mesh->elements->appendData(3, { "TRANSLATION_MOTION", "TRANSLATION_MOTION_X", "TRANSLATION_MOTION_Y", "TRANSLATION_MOTION_Z" });
	}

	if (info::ecf->output.results_selection.gradient || info::ecf->heat_transfer_3d.diffusion_split) {
		_gradient = info::mesh->elements->appendData(3, { "GRADIENT", "GRADIENT_X", "GRADIENT_Y", "GRADIENT_Z" });
	}

	if (info::ecf->output.results_selection.flux) {
		_flux = info::mesh->elements->appendData(3, { "FLUX", "FLUX_X", "FLUX_Y", "FLUX_Z" });
	}

	_boundaries.resize(info::mesh->boundaryRegions.size());
}

HeatTransfer3DController::~HeatTransfer3DController()
{
	delete _kernel;
}

void HeatTransfer3DController::initData()
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

	updateERegions(info::ecf->heat_transfer_3d.initial_temperature, _ntemperature.data->datatarray(), 3, cbegin, tbegin, time);
	updateERegions(_configuration.heat_source, _nheat.data->datatarray(), 3, cbegin, tbegin, time);
	updateERegions(_configuration.translation_motions, _nmotion.data->datatarray(), 3, cbegin, tbegin, time);

	if (info::ecf->heat_transfer_2d.init_temp_respect_bc) {
		initDirichletData(_ntemperature.data->datatarray());
	}
	averageNodeInitilization(_ntemperature.data->datatarray(), _temperature->data);

	if (_motion != NULL) {
		nodeValuesToElements(3, _nmotion.data->datatarray(), _motion->data);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 2) { // TODO: implement edge processing

			auto &distribution = region->procNodes->datatarray().distribution();

			_boundaries[r].coordinate.data = new serializededata<esint, double>(3, distribution);
			_boundaries[r].temperature.data = new serializededata<esint, double>(1, distribution);

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto c = _boundaries[r].coordinate.data->begin(t);
				auto temp = _boundaries[r].temperature.data->begin(t);
				for (auto n = region->procNodes->datatarray().begin(t); n != region->procNodes->datatarray().end(t); ++n, ++c, ++temp) {
					c->at(0) = info::mesh->nodes->coordinates->datatarray()[*n].x;
					c->at(1) = info::mesh->nodes->coordinates->datatarray()[*n].y;
					c->at(2) = info::mesh->nodes->coordinates->datatarray()[*n].z;
					temp->at(0) = _temperature->data[*n];
				}
			}

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();
			tbegin = _boundaries[r].temperature.data->datatarray().begin();

			auto flow = _configuration.heat_flow.find(region->name);
			if (flow != _configuration.heat_flow.end()) {
				updateBRegions(flow->second, _boundaries[r].heatFlow, distribution, 3, cbegin, tbegin, time);
				_boundaries[r].regionArea = info::mesh->boundaryRegions[r]->area;
			}
			auto flux = _configuration.heat_flux.find(region->name);
			if (flux != _configuration.heat_flux.end()) {
				updateBRegions(flux->second, _boundaries[r].heatFlux, distribution, 3, cbegin, tbegin, time);
			}

			auto radiation = _configuration.diffuse_radiation.find(region->name);
			if (radiation != _configuration.diffuse_radiation.end()) {
				updateBRegions(radiation->second.emissivity, _boundaries[r].emissivity, distribution, 3, cbegin, tbegin, time);
				updateBRegions(radiation->second.external_temperature, _boundaries[r].externalTemperature, distribution, 3, cbegin, tbegin, time);
			}

			auto convection = _configuration.convection.find(region->name);
			if (convection != _configuration.convection.end()) {
				if (_boundaries[r].externalTemperature.data == NULL) {
					updateBRegions(convection->second.external_temperature, _boundaries[r].externalTemperature, distribution, 3, cbegin, tbegin, time);
				}
				_boundaries[r].htc.data = new serializededata<esint, double>(1, distribution);

				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
						_boundaries[r].htc.data->datatarray()[i] = _kernel->convectionHTC(convection->second, 3, cbegin + i * 3, time, *(tbegin + i));
					}
				}
			}
		}
	}
}

void HeatTransfer3DController::nextTime()
{
	if (time::isInitial()) {
		return;
	}

	parametersChanged();
}

void HeatTransfer3DController::parametersChanged()
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

	updateERegions(_configuration.heat_source, _nheat.data->datatarray(), 3, cbegin, tbegin, time);
	updateERegions(_configuration.translation_motions, _nmotion.data->datatarray(), 3, cbegin, tbegin, time);

	if (_motion != NULL) {
		nodeValuesToElements(3, _nmotion.data->datatarray(), _motion->data);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (region->dimension == 2) {

			auto &distribution = region->procNodes->datatarray().distribution();

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				auto temp = _boundaries[r].temperature.data->begin(t);
				for (auto n = region->procNodes->datatarray().begin(t); n != region->procNodes->datatarray().end(t); ++n, ++temp) {
					temp->at(0) = _temperature->data[*n];
				}
			}

			cbegin = _boundaries[r].coordinate.data->datatarray().begin();
			tbegin = _boundaries[r].temperature.data->datatarray().begin();

			auto flow = _configuration.heat_flow.find(region->name);
			if (flow != _configuration.heat_flow.end()) {
				updateBRegions(flow->second, _boundaries[r].heatFlow, distribution, 3, cbegin, tbegin, time);
			}
			auto flux = _configuration.heat_flux.find(region->name);
			if (flux != _configuration.heat_flux.end()) {
				updateBRegions(flux->second, _boundaries[r].heatFlux, distribution, 3, cbegin, tbegin, time);
			}

			auto radiation = _configuration.diffuse_radiation.find(region->name);
			if (radiation != _configuration.diffuse_radiation.end()) {
				updateBRegions(radiation->second.emissivity, _boundaries[r].emissivity, distribution, 3, cbegin, tbegin, time);
				updateBRegions(radiation->second.external_temperature, _boundaries[r].externalTemperature, distribution, 3, cbegin, tbegin, time);
			}

			auto convection = _configuration.convection.find(region->name);
			if (convection != _configuration.convection.end()) {
				if (_boundaries[r].externalTemperature.data == NULL) {
					updateBRegions(convection->second.external_temperature, _boundaries[r].externalTemperature, distribution, 3, cbegin, tbegin, time);
				}
				_boundaries[r].htc.data = new serializededata<esint, double>(1, distribution);

				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					for (size_t i = distribution[t]; i < distribution[t + 1]; ++i) {
						_boundaries[r].htc.data->datatarray()[i] = _kernel->convectionHTC(convection->second, 3, cbegin + i * 3, time, *(tbegin + i));
					}
				}
			}
		}
	}
}

void HeatTransfer3DController::processElements(Matrices matrices, const SolverParameters &parameters, InstanceFiller &filler)
{
	auto enodes = info::mesh->elements->procNodes->cbegin() + filler.begin;
	HeatTransfer3DKernel::ElementIterator iterator;

	size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin();
	iterator.temperature = _ntemperature.data->datatarray().begin() + noffset;
	iterator.coordinates = _ncoordinate.data->datatarray().begin() + noffset * 3;
	iterator.motion      = _nmotion.data->datatarray().begin() + noffset * 3;
	iterator.heat        = _nheat.data->datatarray().begin() + noffset;

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

//	if (BEM) {
//		_instance->K[domain].rows = _mesh->domainsSurface->cdistribution[domain + 1] - _mesh->domainsSurface->cdistribution[domain];
//		_instance->K[domain].cols = _instance->K[domain].rows;
//		_instance->K[domain].nnz  = _instance->K[domain].rows * _instance->K[domain].cols;
//		_instance->K[domain].type = 'G';
//		_instance->K[domain].dense_values.resize(_instance->K[domain].nnz);
//		_instance->K[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
//
//		const MaterialConfiguration* material = _mesh->materials[_mesh->procNodes->material->datatarray()[_mesh->procNodes->elementsDistribution[domain]]];
//
//		bem4i::getLaplaceSteklovPoincare(
//				_instance->K[domain].dense_values.data(),
//				_instance->K[domain].rows,
//				reinterpret_cast<double*>(_mesh->domainsSurface->coordinates->datatarray().data() + _mesh->domainsSurface->cdistribution[domain]),
//				_mesh->domainsSurface->tdistribution[domain + 1] - _mesh->domainsSurface->tdistribution[domain],
//				_mesh->domainsSurface->triangles->datatarray().data() + 3 * _mesh->domainsSurface->tdistribution[domain],
//				material->thermal_conductivity.values.get(0, 0).evaluator->evaluate(Point(), _step->currentTime, 0),
//				1,
//				4, 4,
//				_BEM
//	}
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

		size_t noffset = enodes->begin() - info::mesh->elements->procNodes->datatarray().begin(t);
		iterator.temperature = _ntemperature.data->datatarray().begin(t);
		iterator.coordinates = _ncoordinate.data->datatarray().begin(t);
		iterator.motion      = _nmotion.data->datatarray().begin(t);
		iterator.heat        = _nheat.data->datatarray().begin(t);

		if (info::mesh->hasPhaseChange()) {
			iterator.phase = _phaseChange->data.data() + info::mesh->elements->distribution[t];
			iterator.latentHeat = _latentHeat->data.data() + info::mesh->elements->distribution[t];
		}

		if (info::ecf->output.results_selection.gradient) {
			iterator.gradient = _gradient->data.data() + info::mesh->elements->distribution[t] * 3;
		}
		if (info::ecf->output.results_selection.flux) {
			iterator.flux = _flux->data.data() + info::mesh->elements->distribution[t] * 3;
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







