
#include "material.h"

#include "../../configuration.hpp"

espreso::MaterialConfiguration::MaterialConfiguration()
{
	physical_model = PHYSICAL_MODEL::THERMAL;
	REGISTER(physical_model, ECFMetaData()
			.setdescription({ "Select used physical model." })
			.setdatatype({ ECFDataType::OPTION })
			.allowonly([] () { return false; }) // physical model is never printed -> model is set by physics
			.addoption(ECFOption().setname("THERMAL").setdescription("Model used by ADVECTION DIFFUSION."))
			.addoption(ECFOption().setname("LINEAR_ELASTIC").setdescription("One of models used by STRUCTURAL MECHANICS.")));

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

	REGISTER(thermal_properties, ECFMetaData()
            .setdescription({ "Thermal properties" })
			.allowonly([&] () {  return physical_model & PHYSICAL_MODEL::THERMAL; }));
	REGISTER(linear_elastic_properties, ECFMetaData()
            .setdescription({ "Linear elastic properties" })
			.allowonly([&] () {  return physical_model & PHYSICAL_MODEL::LINEAR_ELASTIC; }));
}

