
#include "../instance.h"
#include "../step.h"

#include "../constraints/constraints.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/evaluator/evaluator.h"

#include "../../config/ecf/physics/heattransfer.h"
#include "../assembler/heattransfer.h"

#include "../../mesh/mesh.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/nodestore.h"
#include "../../mesh/store/boundaryregionstore.h"
#include "../../mesh/store/elementsregionstore.h"
#include "../../mesh/store/surfacestore.h"

#include "../../solver/generic/SparseMatrix.h"

using namespace espreso;

HeatTransfer::HeatTransfer(const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration)
: Physics("", NULL, NULL, NULL, &configuration, 1), // skipped because Physics is inherited virtually
  _configuration(configuration), _propertiesConfiguration(propertiesConfiguration),
  _temperature(NULL), _gradient(NULL), _flux(NULL), _phaseChange(NULL), _latentHeat(NULL)
{
	for (eslocal d = 0; d < _mesh->elements->ndomains; d++) {
		if (_BEMDomain[d]) {
			_instance->domainDOFCount[d] = _mesh->domainsSurface->cdistribution[d + 1] - _mesh->domainsSurface->cdistribution[d];
		} else {
			const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[d];
			_instance->domainDOFCount[d] = intervals.back().DOFOffset + intervals.back().end - intervals.back().begin;
		}
	}

	_constraints = new Constraints(
			*_instance,
			*_mesh,
			configuration.load_steps_settings.at(1).temperature,
			configuration.load_steps_settings.at(1).feti.redundant_lagrange,
			configuration.load_steps_settings.at(1).feti.scaling);

	if (_hasBEM) {
		_temperature = _mesh->nodes->appendData(1, { "TEMPERATURE" });
	} else {
		_temperature = _mesh->nodes->appendData(1, { "TEMPERATURE" }, _instance->primalSolution);
	}

	bool phaseChange = false;
	for (size_t m = 0; m < _mesh->materials.size(); m++) {
		phaseChange |= _mesh->materials[m]->phase_change;
	}

	if (phaseChange) {
		_phaseChange = _mesh->nodes->appendData(1, { "PHASE" });
		_latentHeat = _mesh->nodes->appendData(1, { "LATENT_HEAT" });
	}
}

void HeatTransfer::setDirichlet()
{
	if (_step->step) {
		_constraints->update(
				_configuration.load_steps_settings.at(_step->step + 1).temperature,
				_configuration.load_steps_settings.at(_step->step + 1).feti.redundant_lagrange,
				_configuration.load_steps_settings.at(_step->step + 1).feti.scaling);
	}
}

MatrixType HeatTransfer::getMatrixType(size_t domain) const
{
	if (_step->tangentMatrixCorrection) {
		return MatrixType::REAL_UNSYMMETRIC;
	}

	if (_configuration.load_steps_settings.at(_step->step + 1).translation_motions.size()) {
		for (auto it = _configuration.load_steps_settings.at(_step->step + 1).translation_motions.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).translation_motions.end(); ++it) {
			ElementsRegionStore *region = _mesh->eregion(it->first);
			for (eslocal i = _mesh->elements->eintervalsDistribution[domain]; i < _mesh->elements->eintervalsDistribution[domain + 1]; i++) {
				if (region->eintervals[i].begin != region->eintervals[i].end) {
					return MatrixType::REAL_UNSYMMETRIC;
				}
			}
		}
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

void HeatTransfer::prepare()
{

}

void HeatTransfer::analyticRegularization(size_t domain, bool ortogonalCluster)
{
	if (_instance->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC";
	}

	if (
			(_configuration.load_steps_settings.at(_step->step + 1).feti.conjugate_projector != FETI_CONJ_PROJECTOR::CONJ_R &&
			_configuration.load_steps_settings.at(_step->step + 1).feti.conjugate_projector != FETI_CONJ_PROJECTOR::CONJ_K) ||
			_configuration.load_steps_settings.at(_step->step + 1).type != LoadStepConfiguration::TYPE::TRANSIENT) {

		if (_configuration.load_steps_settings.at(_step->step + 1).convection.size()) {
			for (auto it = _configuration.load_steps_settings.at(_step->step + 1).convection.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).convection.end(); ++it) {
				BoundaryRegionStore *region = _mesh->bregion(it->first);
				if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
					return;
				}
			}
		}

		if (_configuration.load_steps_settings.at(_step->step + 1).diffuse_radiation.size()) {
			for (auto it = _configuration.load_steps_settings.at(_step->step + 1).diffuse_radiation.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).diffuse_radiation.end(); ++it) {
				BoundaryRegionStore *region = _mesh->bregion(it->first);
				if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
					return;
				}
			}
		}
	}

	if (_configuration.load_steps_settings.at(_step->step + 1).feti.conjugate_projector == FETI_CONJ_PROJECTOR::CONJ_K) {
		_instance->N1[domain].rows = _instance->K[domain].rows;
		_instance->N1[domain].cols = 1;
		_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
		_instance->N1[domain].type = 'G';

		std::vector<double> diagonal = _instance->K[domain].getDiagonal();
		_instance->N1[domain].dense_values.insert(_instance->N1[domain].dense_values.end(), diagonal.begin(), diagonal.end());

		_instance->RegMat[domain].rows = _instance->K[domain].rows;
		_instance->RegMat[domain].cols = _instance->K[domain].cols;
		_instance->RegMat[domain].nnz  = 1;
		_instance->RegMat[domain].type = _instance->K[domain].type;

		_instance->RegMat[domain].I_row_indices.push_back(1);
		_instance->RegMat[domain].J_col_indices.push_back(1);
		_instance->RegMat[domain].V_values.push_back(_instance->K[domain].getDiagonalMaximum());
		_instance->RegMat[domain].ConvertToCSR(1);
	} else {
		double value;
		if (ortogonalCluster) {
			size_t nSum = 0;
			for (size_t d = 0; d < _instance->domains; d++) {
				if (_mesh->elements->clusters[d] == _mesh->elements->clusters[domain]) {
					nSum += _instance->K[d].rows;
				}
			}
			value = 1 / sqrt(nSum);
		} else {
			value = 1 / sqrt(_instance->K[domain].rows);
		}

		_instance->N1[domain].rows = _instance->K[domain].rows;
		_instance->N1[domain].cols = 1;
		_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
		_instance->N1[domain].type = 'G';

		_instance->N1[domain].dense_values.resize(_instance->N1[domain].nnz, value);

		_instance->RegMat[domain].rows = _instance->K[domain].rows;
		_instance->RegMat[domain].cols = _instance->K[domain].cols;
		_instance->RegMat[domain].nnz  = 1;
		_instance->RegMat[domain].type = _instance->K[domain].type;

		_instance->RegMat[domain].I_row_indices.push_back(1);
		_instance->RegMat[domain].J_col_indices.push_back(1);
		_instance->RegMat[domain].V_values.push_back(_instance->K[domain].getDiagonalMaximum());
		_instance->RegMat[domain].ConvertToCSR(1);
	}
}

void HeatTransfer::computeInitialTemperature(std::vector<std::vector<double> > &data)
{
	auto itemp = _configuration.initial_temperature.find("ALL_ELEMENTS");
	if (_configuration.initial_temperature.size() > 1) {
		ESINFO(GLOBAL_ERROR) << "Not implemented format of initial temperature.";
	}
	if (itemp == _configuration.initial_temperature.end() && _configuration.initial_temperature.size() != 0) {
		ESINFO(GLOBAL_ERROR) << "Not implemented format of initial temperature.";
	}
	if (itemp != _configuration.initial_temperature.end()) {
		if (itemp->second.evaluator->isCoordinateDependent() || itemp->second.evaluator->isTemperatureDependent()) {
			ESINFO(GLOBAL_ERROR) << "Not implemented format of initial temperature.";
		}
	}

	data.resize(_mesh->elements->ndomains);

	#pragma omp parallel for
	for (eslocal d = 0; d < _mesh->elements->ndomains; d++) {
		const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[d];
		data[d].resize(intervals.back().DOFOffset + intervals.back().end - intervals.back().begin, 273.15 + 20);
		if (itemp != _configuration.initial_temperature.end()) {
			itemp->second.evaluator->evalVector(data[d].size(), 0, NULL, NULL, _step->currentTime, data[d].data());
		}
	}

	if (_configuration.init_temp_respect_bc) {
		for (auto it = _configuration.load_steps_settings.at(1).temperature.regions.begin(); it != _configuration.load_steps_settings.at(1).temperature.regions.end(); ++it) {
			BoundaryRegionStore *region = _mesh->bregion(it->first);
			#pragma omp parallel for
			for (eslocal d = 0; d < _mesh->elements->ndomains; d++) {
				const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[d];
				for (size_t i = 0; i < intervals.size(); ++i) {
					for (eslocal n = region->nintervals[intervals[i].pindex].begin; n < region->nintervals[intervals[i].pindex].end; ++n) {
						eslocal dof = intervals[i].DOFOffset + region->nodes->datatarray()[n] - intervals[i].begin;
						data[d][dof] = it->second.evaluator->evaluate(_mesh->nodes->coordinates->datatarray()[region->nodes->datatarray()[n]], 0, 0);
					}
				}
			}
		}
	}

	if (_configuration.diffusion_split) {
		processSolution();
	}
}

void HeatTransfer::preprocessData()
{
	computeInitialTemperature(_instance->primalSolution);
}

double HeatTransfer::sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction) const
{
	switch (restriction) {
	case SumRestriction::NONE:
		return _mesh->sumSquares(data, _mesh->bregion("ALL_NODES"));
	case SumRestriction::DIRICHLET: {
		double sum = 0;
		for (auto it = _configuration.load_steps_settings.at(_step->step + 1).temperature.regions.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).temperature.regions.end(); ++it) {
			sum += _mesh->sumSquares(data, _mesh->bregion(it->first));
		}
		return sum;
	}
	case SumRestriction::NON_DIRICHLET:
		return sumSquares(data, SumRestriction::NONE) - sumSquares(data, SumRestriction::DIRICHLET);
	}
	return 0;
}

void HeatTransfer::convectionMatParameters(
		const ConvectionConfiguration *convection, const Point &p, double temp, double T_EXT,
		double &rho, double &dynamic_viscosity, double &dynamic_viscosity_T, double &heat_capacity, double &thermal_conductivity) const
{
	double  gas_constant;

	switch (convection->fluid) {
	case ConvectionConfiguration::FLUID::AIR:{
		gas_constant = 286.9;
		rho = convection->absolute_pressure.evaluator->evaluate(p, T_EXT, _step->currentTime) / (gas_constant * T_EXT);

		if ((T_EXT >=200) && (T_EXT <= 1600)){
			heat_capacity = 1047.63657-0.372589265*T_EXT+9.45304214E-4*pow(T_EXT,2.0)-6.02409443E-7*pow(T_EXT,3.0)+1.2858961E-10*pow(T_EXT,4.0);
		}else if (T_EXT < 200){
			heat_capacity = 1047.63657-0.372589265*200.0+9.45304214E-4*pow(200.0,2.0)-6.02409443E-7*pow(200.0,3.0)+1.2858961E-10*pow(200.0,4.0);
		}else if (T_EXT > 1600){
			heat_capacity = 1047.63657-0.372589265*1600.0+9.45304214E-4*pow(1600.0,2.0)-6.02409443E-7*pow(1600.0,3.0)+1.2858961E-10*pow(1600.0,4.0);
		}

		if ((T_EXT >=200) && (T_EXT <= 1600)){
			thermal_conductivity = -0.00227583562+1.15480022E-4*T_EXT-7.90252856E-8*pow(T_EXT,2.0)+4.11702505E-11*pow(T_EXT,3.0)-7.43864331E-15*pow(T_EXT,4.0);
		}else if (T_EXT < 200){
			thermal_conductivity = -0.00227583562+1.15480022E-4*200.0-7.90252856E-8*pow(200.0,2.0)+4.11702505E-11*pow(200.0,3.0)-7.43864331E-15*pow(200.0,4.0);
		}else if (T_EXT > 1600){
			thermal_conductivity =  -0.00227583562+1.15480022E-4*1600.0-7.90252856E-8*pow(1600.0,2.0)+4.11702505E-11*pow(1600.0,3.0)-7.43864331E-15*pow(1600.0,4.0);
		}

		if ((T_EXT >=200) && (T_EXT <= 1600)){
		dynamic_viscosity = -8.38278E-7 + 8.35717342E-8 * T_EXT - 7.69429583E-11 * pow(T_EXT,2.0) + 4.6437266E-14 * pow(T_EXT,3.0) - 1.06585607E-17 * pow(T_EXT,4.0);
		}else if (T_EXT < 200){
			dynamic_viscosity = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * pow(200.0,2.0) + 4.6437266E-14 * pow(200.0,3.0) - 1.06585607E-17 * pow(200.0,4.0);
		}else if (T_EXT > 1600){
			dynamic_viscosity = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * pow(1600.0,2.0) + 4.6437266E-14 * pow(1600.0,3.0) - 1.06585607E-17 * pow(1600.0,4.0);
		}


		if ((temp >=200) && (temp <= 1600)){
			dynamic_viscosity_T = -8.38278E-7 + 8.35717342E-8 * temp - 7.69429583E-11 * pow(temp,2.0) + 4.6437266E-14 * pow(temp,3.0) - 1.06585607E-17 * pow(temp,4.0);
		}else if (temp < 200){
			dynamic_viscosity_T = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * pow(200.0,2.0) + 4.6437266E-14 * pow(200.0,3.0) - 1.06585607E-17 * pow(200.0,4.0);
		}else if (temp > 1600){
			dynamic_viscosity_T = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * pow(1600.0,2.0) + 4.6437266E-14 * pow(1600.0,3.0) - 1.06585607E-17 * pow(1600.0,4.0);
		}

	}break;
	case ConvectionConfiguration::FLUID::WATER:{

		if ((T_EXT >=273) && (T_EXT <= 283)){
			rho = 972.7584 + 0.2084 *T_EXT - 4.0E-4 * pow(T_EXT,2.0);
		}else if((T_EXT >283) && (T_EXT <= 373)){
			rho = 345.28 + 5.749816 * T_EXT - 0.0157244 * pow(T_EXT,2.0) + 1.264375E-5 * pow(T_EXT,3.0);
		}else if (T_EXT < 273){
			rho = 972.7584 + 0.2084 *273.0 - 4.0E-4 * pow(273.0,2.0);
		}else if (T_EXT > 373){
			rho = 345.28 + 5.749816 * 373.0 - 0.0157244 * pow(373.0,2.0) + 1.264375E-5 * pow(373.0,3.0);
		}


		if ((T_EXT >= 265) && (T_EXT <= 293)){
			dynamic_viscosity = 5.948859 - 0.08236196 * T_EXT + 4.287142E-4 * pow(T_EXT,2.0) - 9.938045E-7 * pow(T_EXT,3.0) + 8.65316E-10 * pow(T_EXT,4.0);
		}else if((T_EXT >293) && (T_EXT <= 353)){
			dynamic_viscosity = 	0.410191 - 0.004753985 * T_EXT + 2.079795E-5 * pow(T_EXT,2.0) - 4.061698E-8 *  pow(T_EXT,3.0) + 2.983925E-11 * pow(T_EXT,4.0);
		}else if((T_EXT >353) && (T_EXT <= 423)){
			dynamic_viscosity = 0.03625638 - 3.265463E-4 * T_EXT + 1.127139E-6 * pow(T_EXT,2.0) - 1.75363E-9 * pow(T_EXT,3.0) + 1.033976E-12 * pow(T_EXT,4.0);
		}else if (T_EXT < 265){
			dynamic_viscosity = 5.948859 - 0.08236196 * 265.0 + 4.287142E-4 * pow(265.0,2.0) - 9.938045E-7 * pow(265.0,3.0) + 8.65316E-10 * pow(265.0,4.0);
		}else if (T_EXT > 423){
			dynamic_viscosity = 0.03625638 - 3.265463E-4 * 423.0 + 1.127139E-6 * pow(423.0,2.0) - 1.75363E-9 * pow(423.0,3.0) + 1.033976E-12 * pow(423.0,4.0);
		}


		if ((T_EXT >= 275) && (T_EXT <= 370)){
			thermal_conductivity = -0.9003748 + 0.008387698 * T_EXT - 1.118205E-5 * pow(T_EXT,2.0);
		}else if (T_EXT < 275){
			thermal_conductivity = -0.9003748 + 0.008387698 * 275.0 - 1.118205E-5 * pow(275.0,2.0);
		}else if (T_EXT > 370){
			thermal_conductivity = -0.9003748 + 0.008387698 * 370.0 - 1.118205E-5 * pow(370.0,2.0);
		}

		if ((T_EXT >= 293) && (T_EXT <= 373)){
			heat_capacity = 4035.841 + 0.492312 * T_EXT;
		}else if (T_EXT < 293){
			heat_capacity = 4035.841 + 0.492312 * 293.0;
		}else if (T_EXT > 373){
			heat_capacity = 4035.841 + 0.492312 * 373.0;
		}


		if ((temp >= 265) && (temp <= 293)){
			dynamic_viscosity_T = 5.948859 - 0.08236196 * temp + 4.287142E-4 * pow(temp,2.0) - 9.938045E-7 * pow(temp,3.0) + 8.65316E-10 * pow(temp,4.0);
		}else if((temp >293) && (temp <= 353)){
			dynamic_viscosity_T = 	0.410191 - 0.004753985 * temp + 2.079795E-5 * pow(temp,2.0) - 4.061698E-8 *  pow(temp,3.0) + 2.983925E-11 * pow(temp,4.0);
		}else if((temp >353) && (temp <= 423)){
			dynamic_viscosity_T = 0.03625638 - 3.265463E-4 * temp + 1.127139E-6 * pow(temp,2.0) - 1.75363E-9 * pow(temp,3.0) + 1.033976E-12 * pow(temp,4.0);
		}else if (temp < 265){
			dynamic_viscosity_T = 5.948859 - 0.08236196 * 265.0 + 4.287142E-4 * pow(265.0,2.0) - 9.938045E-7 * pow(265.0,3.0) + 8.65316E-10 * pow(265.0,4.0);
		}else if (temp > 423){
			dynamic_viscosity_T = 0.03625638 - 3.265463E-4 * 423.0 + 1.127139E-6 * pow(423.0,2.0) - 1.75363E-9 * pow(423.0,3.0) + 1.033976E-12 * pow(423.0,4.0);
		}

	}break;
	case ConvectionConfiguration::FLUID::ENGINE_OIL:{


		if ((T_EXT >=273) && (T_EXT <= 433)){
			rho = 1068.70404 - 0.6393421 * T_EXT + 7.34307359E-5 * pow(T_EXT,2.0);
		}else if (T_EXT < 273){
			rho = 1068.70404 - 0.6393421 * 273.0 + 7.34307359E-5 * pow(273.0,2.0);
		}else if (T_EXT > 433){
			rho = 1068.70404 - 0.6393421 * 433.0 + 7.34307359E-5 * pow(433.0,2.0);
		}

		if ((T_EXT >= 273) && (T_EXT <= 433)){
			thermal_conductivity = 0.192223542 - 2.0637987E-4 * T_EXT + 1.54220779E-7 * pow(T_EXT,2.0);
		}else if (T_EXT < 273){
			thermal_conductivity = 0.192223542 - 2.0637987E-4 * 273.0 + 1.54220779E-7 * pow(273.0,2.0);
		}else if (T_EXT > 433){
			thermal_conductivity =0.192223542 - 2.0637987E-4 * 433.0 + 1.54220779E-7 * pow(433.0,2.0);
		}


		if ((T_EXT >= 273) && (T_EXT <= 433)){
			heat_capacity = 761.405625 + 3.47685606 * T_EXT + 0.00115530303 * pow(T_EXT,2.0);
		}else if (T_EXT < 273){
			heat_capacity = 761.405625 + 3.47685606 * 273.0 + 0.00115530303 * pow(273.0,2.0);
		}else if (T_EXT > 433){
			heat_capacity = 761.405625 + 3.47685606 * 433.0 + 0.00115530303 * pow(433.0,2.0);
		}


		if ((T_EXT >= 273) && (T_EXT <= 353)){
			dynamic_viscosity = 42669.28688622 - 741.1718801282 * T_EXT + 5.360521287088 * pow(T_EXT,2.0) - 0.02066027676164 * pow(T_EXT,3.0) + 4.47491538052E-5 * pow(T_EXT,4.0) - 5.164053479202E-8 * pow(T_EXT,5.0) + 2.48033770504E-11 * pow(T_EXT,6.0);
		}else if ((T_EXT > 353) && (T_EXT <= 433 )){
			dynamic_viscosity = 4.94593941 - 0.0351869631 * T_EXT + 8.37935977E-5 * pow(T_EXT,2.0) - 6.67125E-8 * pow(T_EXT,3.0);

		}else if (T_EXT < 273){
			dynamic_viscosity = 42669.28688622 - 741.1718801282 * 273.0 + 5.360521287088 * pow(273.0,2.0) - 0.02066027676164 * pow(273.0,3.0) + 4.47491538052E-5 * pow(273.0,4.0) - 5.164053479202E-8 * pow(273.0,5.0) + 2.48033770504E-11 * pow(273.0,6.0);
		}else if (T_EXT > 433){
			dynamic_viscosity = 4.94593941 - 0.0351869631 * 433.0 + 8.37935977E-5 * pow(433.0,2.0) - 6.67125E-8 * pow(433.0,3.0);
		}

		if ((temp >= 273) && (temp <= 353)){
			dynamic_viscosity_T = 42669.28688622 - 741.1718801282 * temp + 5.360521287088 * pow(temp,2.0) - 0.02066027676164 * pow(temp,3.0) + 4.47491538052E-5 * pow(temp,4.0) - 5.164053479202E-8 * pow(temp,5.0) + 2.48033770504E-11 * pow(temp,6.0);
		}else if ((temp > 353) && (temp <= 433 )){
			dynamic_viscosity_T = 4.94593941 - 0.0351869631 * temp + 8.37935977E-5 * pow(temp,2.0) - 6.67125E-8 * pow(temp,3.0);

		}else if (temp < 273){
			dynamic_viscosity_T = 42669.28688622 - 741.1718801282 * 273.0 + 5.360521287088 * pow(273.0,2.0) - 0.02066027676164 * pow(273.0,3.0) + 4.47491538052E-5 * pow(273.0,4.0) - 5.164053479202E-8 * pow(273.0,5.0) + 2.48033770504E-11 * pow(273.0,6.0);
		}else if (temp > 433){
			dynamic_viscosity_T = 4.94593941 - 0.0351869631 * 433.0 + 8.37935977E-5 * pow(433.0,2.0) - 6.67125E-8 * pow(433.0,3.0);
		}


	}break;
	case ConvectionConfiguration::FLUID::TRANSFORMER_OIL:{

		if ((T_EXT >=223) && (T_EXT <= 373)){
			rho = 1055.04607 - 0.581753034 * T_EXT - 6.40531689E-5 * pow(T_EXT,2.0);
		}else if (T_EXT < 223){
			rho =  1055.04607 - 0.581753034 * 223.0 - 6.40531689E-5 * pow(223.0,2.0);
		}else if (T_EXT > 373){
			rho = 1055.04607 - 0.581753034 * 373.0 - 6.40531689E-5 * pow(373.0,2.0);
		}

		if ((T_EXT >= 273) && (T_EXT <= 433)){
			thermal_conductivity = 0.134299084 - 8.04973822E-5 * T_EXT;
		}else if (T_EXT < 273){
			thermal_conductivity = 0.134299084 - 8.04973822E-5 * 223.0;
		}else if (T_EXT > 433){
			thermal_conductivity = 0.134299084 - 8.04973822E-5 * 373.0;
		}


		if ((T_EXT >= 223) && (T_EXT <= 293)){
			heat_capacity = -117056.38 + 1816.76208 * T_EXT - 10.305786 * pow(T_EXT,2.0) + 0.0256691919 * pow(T_EXT,3.0) - 2.36742424E-5 * pow(T_EXT,4.0);
		}else if ((T_EXT > 293) && (T_EXT <= 373 )){
			heat_capacity = -13408.1491 + 123.044152 * T_EXT - 0.335401786 * pow(T_EXT,2.0) + 3.125E-4 * pow(T_EXT,3.0);
		}else if (T_EXT < 223){
			heat_capacity = -117056.38 + 1816.76208 * 223.0 - 10.305786 * pow(223.0,2.0) + 0.0256691919 * pow(223.0,3.0) - 2.36742424E-5 * pow(223.0,4.0);
		}else if (T_EXT > 373){
			heat_capacity = -13408.1491 + 123.044152 * 373.0 - 0.335401786 * pow(373.0,2.0) + 3.125E-4 * pow(373.0,3.0);
		}


		if ((T_EXT >= 243) && (T_EXT <= 273)){
			dynamic_viscosity = 4492.20229 - 64.7408879 * T_EXT + 0.349900959 * pow(T_EXT,2.0) - 8.40477E-4 * pow(T_EXT,3.0) + 7.57041667E-7 * pow(T_EXT,4.0);
		}else if ((T_EXT > 273) && (T_EXT <= 373 )){
			dynamic_viscosity = 91.4524999 - 1.33227058 * T_EXT + 0.00777680216 * pow(T_EXT,2.0) - 2.27271368E-5 *  pow(T_EXT,3.0) + 3.32419673E-8 * pow(T_EXT,4.0) - 1.94631023E-11 * pow(T_EXT,5.0);
		}else if (T_EXT < 243){
			dynamic_viscosity = 4492.20229 - 64.7408879 * 243.0 + 0.349900959 * pow(243.0,2.0) - 8.40477E-4 * pow(243.0,3.0) + 7.57041667E-7 * pow(243.0,4.0);
		}else if (T_EXT > 373){
			dynamic_viscosity = 91.4524999 - 1.33227058 * 373.0 + 0.00777680216 * pow(373.0,2.0) - 2.27271368E-5 *  pow(373.0,3.0) + 3.32419673E-8 * pow(373.0,4.0) - 1.94631023E-11 * pow(373.0,5.0);

		}

		if ((temp >= 243) && (temp <= 273)){
			dynamic_viscosity_T = 4492.20229 - 64.7408879 * temp + 0.349900959 * pow(temp,2.0) - 8.40477E-4 * pow(temp,3.0) + 7.57041667E-7 * pow(temp,4.0);
		}else if ((temp > 273) && (temp <= 373 )){
			dynamic_viscosity_T = 91.4524999 - 1.33227058 * temp + 0.00777680216 * pow(temp,2.0) - 2.27271368E-5 *  pow(temp,3.0) + 3.32419673E-8 * pow(temp,4.0) - 1.94631023E-11 * pow(temp,5.0);
		}else if (temp < 243){
			dynamic_viscosity_T = 4492.20229 - 64.7408879 * 243.0 + 0.349900959 * pow(243.0,2.0) - 8.40477E-4 * pow(243.0,3.0) + 7.57041667E-7 * pow(243.0,4.0);
		}else if (temp > 373){
			dynamic_viscosity_T = 91.4524999 - 1.33227058 * 373.0 + 0.00777680216 * pow(373.0,2.0) - 2.27271368E-5 *  pow(373.0,3.0) + 3.32419673E-8 * pow(373.0,4.0) - 1.94631023E-11 * pow(373.0,5.0);

		}

	}break;
	default:
		ESINFO(ERROR) << "Invalid convection fluid type.";
	}
}


double HeatTransfer::computeHTC(const ConvectionConfiguration *convection, const Point &p, double temp) const
{
	double htc = 0;
	switch (convection->type) {
	case ConvectionConfiguration::TYPE::USER:
		return convection->heat_transfer_coefficient.evaluator->evaluate(p, temp, _step->currentTime);

	case ConvectionConfiguration::TYPE::EXTERNAL_NATURAL: {
		double T_AVG, g, rho, dynamic_viscosity, heat_capacity, thermal_conductivity, dynamic_viscosity_T, text;

		text = convection->external_temperature.evaluator->evaluate(p, temp, _step->currentTime);
		T_AVG = (text + temp) / 2;
		g = 9.81;

		convectionMatParameters(convection, p, temp, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity);

		switch (convection->variant) {
		case ConvectionConfiguration::VARIANT::INCLINED_WALL: {
			double wallHeaight = convection->wall_height.evaluator->evaluate(p, temp, _step->currentTime);
			double RaL = pow(rho,2) * g * (1 / T_AVG) * heat_capacity * std::fabs(temp - text  ) * pow(wallHeaight, 3.0) / ( thermal_conductivity * dynamic_viscosity);
			double tilt_angle = convection->tilt_angle.evaluator->evaluate(p, temp, _step->currentTime) * M_PI / 180.0;
			if (RaL <= 1e9) {
				htc = (thermal_conductivity / wallHeaight) * (0.68 + (0.67 * cos(tilt_angle) * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermal_conductivity / wallHeaight) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
			}

		} break;
		case ConvectionConfiguration::VARIANT::VERTICAL_WALL: {
			double wallHeaight = convection->wall_height.evaluator->evaluate(p, temp, _step->currentTime);
			double RaL = pow(rho,2) * g * (1/T_AVG) * heat_capacity * std::fabs(temp - text) * pow(wallHeaight, 3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (RaL <= 1e9) {
				htc = (thermal_conductivity / wallHeaight) * (0.68 + (0.67 * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermal_conductivity / wallHeaight) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
			}

		} break;
		case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_UP: {
			double lenght = convection->length.evaluator->evaluate(p, temp, _step->currentTime);
			double RaL = pow(rho,2) * g * (1/T_AVG) * heat_capacity * std::fabs(temp - text) * pow(lenght, 3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (temp > text) {
				if (RaL <= 1e7) {
					htc = thermal_conductivity / lenght * 0.54 * pow(RaL,0.25);
				} else {
					htc = thermal_conductivity / lenght * 0.15 * pow(RaL,1.0/3.0);
				}
			} else {
				htc = thermal_conductivity / lenght * 0.27 * pow(RaL,0.25);
			}

		} break;
		case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_DOWN: {
			double lenght = convection->length.evaluator->evaluate(p, temp, _step->currentTime);
			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs(temp - text) *pow(lenght, 3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (temp <= text) {
				if (RaL <= 1e7) {
					htc = thermal_conductivity / lenght * 0.54 * pow(RaL,0.25);
				} else {
					htc = thermal_conductivity / lenght * 0.15 * pow(RaL,1.0/3.0);
				}
			} else {
				htc = thermal_conductivity / lenght * 0.27 * pow(RaL,0.25);
			}
		} break;
		case ConvectionConfiguration::VARIANT::HORIZONTAL_CYLINDER:{
			double diameter = convection->diameter.evaluator->evaluate(p, temp, _step->currentTime);
			double RaD = pow(rho,2) * g * (1/T_AVG) * heat_capacity * std::fabs(temp - text) * pow(diameter, 3.0)/ ( thermal_conductivity * dynamic_viscosity);
			double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;

			if ( RaD > 10e12 ){
				// warning!!!!
				ESINFO(ERROR) << "Validated only for RaD <= 10e12 ";
			}
			htc = thermal_conductivity / diameter * pow( 0.6 + ( 0.387*pow(RaD,1.0/6.0)/ pow( 1 + pow( 0.559/Pr, 9.0/16.0), 8.0/27.0) ) ,2.0);

		} break;
		case ConvectionConfiguration::VARIANT::SPHERE: {
			double diameter = convection->diameter.evaluator->evaluate(p, temp, _step->currentTime);
			double RaD = pow(rho,2) * g * (1/T_AVG) * heat_capacity * std::fabs( temp - text) * pow(diameter, 3.0) / ( thermal_conductivity * dynamic_viscosity);
			double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;

			if ( RaD > 10e11 || Pr < 0.7 ){
				// warning!!!!
				ESINFO(ERROR) << "Validated only for RaD <= 10e11 and Pr >= 0.7 ";
			}
			htc = thermal_conductivity / diameter * pow(2.0 + ( 0.589*pow(RaD,0.25)/ pow( 1 + pow( 0.469/Pr, 9.0/16.0), 4.0/9.0) ) ,2.0);
		} break;
		default:
			ESINFO(ERROR) << "Invalid convection variant for EXTERNAL_NATURAL.";
		}
	} break;

	case ConvectionConfiguration::TYPE::INTERNAL_NATURAL:{

		double T_AVG, g, rho, dynamic_viscosity, heat_capacity, thermal_conductivity,dynamic_viscosity_T, text;

		text = convection->external_temperature.evaluator->evaluate(p, temp, _step->currentTime);
		T_AVG = (text + temp) / 2.0;
		g = 9.81;

		convectionMatParameters(convection, p, temp, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity );

		switch (convection->variant) {
		case ConvectionConfiguration::VARIANT::PARALLEL_PLATES: {
			double wallHeight = convection->wall_height.evaluator->evaluate(p, temp, _step->currentTime);
			double length = convection->length.evaluator->evaluate(p, temp, _step->currentTime);
			double H_L = wallHeight / length;
			double RaL = pow(rho,2) * g * (1/T_AVG) * heat_capacity * std::fabs(temp - text) * pow(length,3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (( RaL < H_L ) && (temp >  text)) {
				htc = thermal_conductivity / wallHeight * ( 1.0 / 24.0 ) * RaL;
			} else {
				double RaL = pow(rho,2) * g * (1/T_AVG) * heat_capacity * std::fabs( temp - text) * pow(length,3.0)/ ( thermal_conductivity * dynamic_viscosity);
				if (RaL <= 1e9) {
					htc = (thermal_conductivity / length) * (0.68 + (0.67 * pow(RaL,0.25)) / (pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
				} else {
					htc = (thermal_conductivity / length) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
				}
			}

		} break;
		case ConvectionConfiguration::VARIANT::CIRCULAR_TUBE: {
			double diameter = convection->diameter.evaluator->evaluate(p, temp, _step->currentTime);
			double wallHeight = convection->wall_height.evaluator->evaluate(p, temp, _step->currentTime);
			double RaD = pow(rho,2) * g * (1/T_AVG) * heat_capacity * std::fabs(temp - text) * pow(diameter, 3.0)/ ( thermal_conductivity * dynamic_viscosity);
			double H_D = wallHeight / diameter;

			if ( RaD < H_D ){
				htc = thermal_conductivity / wallHeight * ( 1.0 / 128.0 ) * RaD;
			}else{

				double RaD = pow(rho,2) * g * (1/T_AVG) * heat_capacity * std::fabs(temp - text) * pow(diameter, 3.0)/ ( thermal_conductivity * dynamic_viscosity);
				if (RaD <= 1e9) {
					htc = (thermal_conductivity / diameter) * (0.68 + (0.67 * pow(RaD,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
				} else {
					htc = (thermal_conductivity / diameter) * pow(0.825 + (0.387 * pow(RaD,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
				}
			}

		} break;
		default:
			ESINFO(ERROR) << "Invalid convection variant for INTERNAL_NATURAL.";
		}
	} break;

	case ConvectionConfiguration::TYPE::EXTERNAL_FORCED: {

			switch (convection->variant) {
			case ConvectionConfiguration::VARIANT::AVERAGE_PLATE: {

				double T_AVG, rho, dynamic_viscosity, heat_capacity, thermal_conductivity,dynamic_viscosity_T, text;
				text = convection->external_temperature.evaluator->evaluate(p, temp, _step->currentTime);
				double length = convection->length.evaluator->evaluate(p, temp, _step->currentTime);

				T_AVG = (text + temp) / 2.0;

				convectionMatParameters(convection, p, temp, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity );

				double Re = rho * convection->fluid_velocity.evaluator->evaluate(p, temp, _step->currentTime) * length / dynamic_viscosity;
				double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;
				if (Re <= 5e5) {
					htc = 2 * (thermal_conductivity / length) * ((0.3387 * pow(Pr, 1.0 / 3.0) * pow(Re, 0.5)) / (pow(1 + pow(0.0468 / Pr, 2.0 / 3.0), 0.25)));
				} else {
					htc = 2 * (thermal_conductivity / length) * pow(Pr, 1.0 / 3.0)	* (0.037 * pow(Re, 0.8) - 871);
				}

			} break;
			default:
				ESINFO(ERROR) << "Invalid convection variant for EXTERNAL_FORCED.";
			}
	} break;

	case ConvectionConfiguration::TYPE::INTERNAL_FORCED:{

			switch (convection->variant) {
			case ConvectionConfiguration::VARIANT::TUBE: {

				double T_EXT, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity;

				T_EXT = convection->external_temperature.evaluator->evaluate(p, temp, _step->currentTime);

				convectionMatParameters(convection, p, temp, T_EXT, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity);

				double Re = rho * convection->fluid_velocity.evaluator->evaluate(p, temp, _step->currentTime) * convection->diameter.evaluator->evaluate(p, temp, _step->currentTime) / dynamic_viscosity;
				double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;
				double n = temp < T_EXT ? 0.3 : 0.4;
				htc = thermal_conductivity / convection->diameter.evaluator->evaluate(p, temp, _step->currentTime);
				if (Re <= 2500) {
					htc *= 3.66;
				} else {
					htc *= 0.027 * pow(Re, .8) * pow(Pr, n)	* pow(dynamic_viscosity / dynamic_viscosity_T, 0.14);
				}
			}break;
			default:
				ESINFO(ERROR) << "Invalid convection variant for EXTERNAL_FORCED.";
			}
	}break;

	default:
		ESINFO(ERROR) << "Invalid convection TYPE.";
	}

	return htc;
}






