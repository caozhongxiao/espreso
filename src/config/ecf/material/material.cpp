
#include "material.h"

#include "../../configuration.hpp"

espreso::MaterialConfiguration::MaterialConfiguration()
: physical_model(static_cast<PHYSICAL_MODEL>(~0)), _allowed_physical_models(static_cast<PHYSICAL_MODEL>(~0))
{
	REGISTER(coordinate_system, ECFMetaData()
			.setdescription({ "Material coordinate system." }));

	REGISTER(physical_model, ECFMetaData()
			.setdescription({ "Select used physical model." })
			.setdatatype({ ECFDataType::ENUM_FLAGS })
			.addoption(ECFOption().setname("THERMAL").setdescription("Model used by HEAT TRANSFER.")
					.allowonly([&] () { return _allowed_physical_models & PHYSICAL_MODEL::THERMAL; }))
			.addoption(ECFOption().setname("LINEAR_ELASTIC").setdescription("One of models used by STRUCTURAL MECHANICS.")
					.allowonly([&] () { return _allowed_physical_models & PHYSICAL_MODEL::LINEAR_ELASTIC; })));

	density.value = heat_capacity.value = "0";
	registerParameter("dens", density, ECFMetaData()
			.setdescription({ "Density" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setunit("kg/m^3")
			.setmaterialvariables());

	registerParameter("CP", heat_capacity, ECFMetaData()
			.setdescription({ "Heat capacity" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setunit("J / (kg * K)")
			.setmaterialvariables());

	REGISTER(thermal_conductivity, ECFMetaData()
			.setdescription({ "Thermal conductivity." })
			.allowonly([&] () {  return physical_model & PHYSICAL_MODEL::THERMAL; }));

	REGISTER(linear_elastic_properties, ECFMetaData()
			.setdescription({ "Linear elastic properties." })
			.allowonly([&] () {  return physical_model & PHYSICAL_MODEL::LINEAR_ELASTIC; }));
}

espreso::MaterialConfiguration::MaterialConfiguration(DIMENSION dimension, PHYSICAL_MODEL allowedPhysicalModels)
: espreso::MaterialConfiguration()
{
	_allowed_physical_models = allowedPhysicalModels;
	coordinate_system.dimension = dimension;
	thermal_conductivity.dimension = dimension;
	linear_elastic_properties.dimension = dimension;

	// Material has special behavior.
	// GUI can edit materials for all physics.
	// Default physics should print only subset of parameters.

	// DROP not allowed parameters.
	physical_model = allowedPhysicalModels;
	for (size_t i = 0; i < parameters.size();) {
		if (!parameters[i]->metadata.isallowed()) {
			dropParameter(parameters[i]);
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
