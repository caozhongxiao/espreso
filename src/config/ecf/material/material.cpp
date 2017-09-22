
#include "material.h"

#include "../../configuration.hpp"

using namespace espreso;
: _phase_change(NULL), _physical_model(NULL)
{
	REGISTER(coordinate_system, ECFMetaData()
			.setdescription({ "Coordinate system" })
			.allowonly([&] () { return _phase_change == NULL || !*_phase_change; }));

	density.value = heat_capacity.value = "0";
	registerParameter("dens", density, ECFMetaData()
			.setdescription({ "Density" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setunit("kg/m^3")
			.setmaterialvariables()
			.allowonly([&] () { return _phase_change == NULL || !*_phase_change; }));

	registerParameter("CP", heat_capacity, ECFMetaData()
			.setdescription({ "Heat capacity" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setunit("J / (kg * K)")
			.setmaterialvariables()
			.allowonly([&] () { return _phase_change == NULL || !*_phase_change; }));

	REGISTER(thermal_conductivity, ECFMetaData()
			.setdescription({ "Thermal conductivity" })
			.allowonly([&] () {
				return (_phase_change == NULL || !*_phase_change) && (*_physical_model & PHYSICAL_MODEL::THERMAL);
			}));

	REGISTER(linear_elastic_properties, ECFMetaData()
			.setdescription({ "Linear elasticity" })
			.allowonly([&] () { return (_phase_change == NULL || !*_phase_change) && (*_physical_model & PHYSICAL_MODEL::LINEAR_ELASTIC); }));
}

MaterialBaseConfiguration::MaterialBaseConfiguration(bool *phase_change, PHYSICAL_MODEL *physicalModel)
: MaterialBaseConfiguration()
{
	_phase_change = phase_change;
	_physical_model = physicalModel;
}

MaterialBaseConfiguration::MaterialBaseConfiguration(bool *phase_change, PHYSICAL_MODEL *physicalModel, DIMENSION *dimension)
: MaterialBaseConfiguration(phase_change, physicalModel)
{
	coordinate_system.dimension = *dimension;
	thermal_conductivity.dimension = *dimension;
	linear_elastic_properties.dimension = *dimension;

	// DROP not allowed parameters.
	for (size_t i = 0; i < parameters.size();) {
		if (!parameters[i]->metadata.isallowed()) {
			dropParameter(parameters[i]);
		} else {
			i++;
		}
	}
}

MaterialConfiguration::MaterialConfiguration()
: MaterialBaseConfiguration(&phase_change, &physical_model),
  dimension(DIMENSION::D3),
  physical_model(static_cast<PHYSICAL_MODEL>(~0)),
  phase_change(false),
  _allowed_physical_models(static_cast<PHYSICAL_MODEL>(~0))
{
	name = "MATERIAL_NAME";
	REGISTER(name, ECFMetaData()
			 .setdescription({ "Name" })
			 .setdatatype( { ECFDataType::STRING } ));

	description = "MATERIAL_DESCRIPTION";
	REGISTER(description, ECFMetaData()
			 .setdescription({ "Description" })
			 .setdatatype( { ECFDataType::STRING } ));

	addSpace();

    physical_model = PHYSICAL_MODEL::THERMAL;
	REGISTER(physical_model, ECFMetaData()
			.setdescription({ "Physical model" })
			.setdatatype({ ECFDataType::ENUM_FLAGS })
			.addoption(ECFOption().setname("THERMAL").setdescription("Model used by HEAT TRANSFER.")
					.allowonly([&] () { return _allowed_physical_models & PHYSICAL_MODEL::THERMAL; }))
			.addoption(ECFOption().setname("LINEAR_ELASTIC").setdescription("One of models used by STRUCTURAL MECHANICS.")
					.allowonly([&] () { return _allowed_physical_models & PHYSICAL_MODEL::LINEAR_ELASTIC; })));

	REGISTER(phase_change, ECFMetaData()
			.setdescription({ "Turn on/off phase change." })
			.setdatatype({ ECFDataType::BOOL }));

	addSeparator();

	moveLastBefore(parameters.front()->name);
	moveLastBefore(parameters.front()->name);
	moveLastBefore(parameters.front()->name);
	moveLastBefore(parameters.front()->name);
	moveLastBefore(parameters.front()->name);
	moveLastBefore(parameters.front()->name);

	smooth_step_order = 1;
	REGISTER(smooth_step_order, ECFMetaData()
			.setdescription({ "Smooth step order." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER })
			.allowonly([&] () { return phase_change; }));

	latent_heat = transition_interval = phase_change_temperature = 0;
	REGISTER(latent_heat, ECFMetaData()
			.setdescription({ "Latent heat." })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return phase_change; }));
	REGISTER(transition_interval, ECFMetaData()
				.setdescription({ "Transition interval." })
				.setdatatype({ ECFDataType::FLOAT })
				.allowonly([&] () { return phase_change; }));
	REGISTER(phase_change_temperature, ECFMetaData()
				.setdescription({ "Phase change temperature." })
				.setdatatype({ ECFDataType::FLOAT })
				.allowonly([&] () { return phase_change; }));

	REGISTER(phases, ECFMetaData()
			.setdescription({ "Phase change settings.", "Settings for a particular phase." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER })
			.setpattern({ "1" })
			.allowonly([&] () { return phase_change; }),
			static_cast<bool*>(NULL), &physical_model, &dimension);

	getParameter(&phases)->getParameter("1");
	getParameter(&phases)->getParameter("2");
}

MaterialConfiguration::MaterialConfiguration(DIMENSION dimension, PHYSICAL_MODEL allowedPhysicalModels)
: MaterialConfiguration()
{
	this->dimension = dimension;
	_allowed_physical_models = allowedPhysicalModels;
	coordinate_system.dimension = dimension;
	thermal_conductivity.dimension = dimension;
	linear_elastic_properties.dimension = dimension;

	// Material has special behavior.
	// GUI can edit materials for all physics.
	// Default physics should print only subset of parameters.

	// DROP not allowed models
	physical_model = allowedPhysicalModels;
	bool phaseChange = phase_change;
	for (size_t i = 0; i < parameters.size();) {
		if (!parameters[i]->metadata.isallowed()) {
			phase_change = !phase_change;
			if (!parameters[i]->metadata.isallowed()) {
				dropParameter(parameters[i]);
			} else {
				i++;
			}
			phase_change = phaseChange;
		} else {
			i++;
		}
	}

	// 2. Set physical_model to the first allowed.
	for (int i = 0; i < 64; i++) {
		if ((1 << i) >= static_cast<int>(allowedPhysicalModels)) {
			physical_model = static_cast<PHYSICAL_MODEL>(1 << i);
			break;
		}
	}
}

MaterialConfiguration& MaterialConfiguration::operator=(const MaterialConfiguration &other)
{
	if (this != &other) {
		name = other.name;
		description = other.description;
		dimension = other.dimension;
		physical_model = other.physical_model;
		_allowed_physical_models = other._allowed_physical_models;

		phase_change = other.phase_change;
		smooth_step_order = other.smooth_step_order;
		latent_heat = other.latent_heat;
		transition_interval = other.transition_interval;
		phase_change_temperature = other.phase_change_temperature;

		*dynamic_cast<MaterialBaseConfiguration*>(this) = dynamic_cast<const MaterialBaseConfiguration&>(other);

		for (auto it = other.phases.begin(); it != other.phases.end(); ++it) {
			getParameter(&phases)->getParameter(std::to_string(it->first));
			phases[it->first] = it->second;
		}
	}
	return *this;
}
