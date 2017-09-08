
#include "material.h"

#include "../../configuration.hpp"

espreso::MaterialConfiguration::MaterialConfiguration()
: _allowed_physical_models(static_cast<PHYSICAL_MODEL>(~0))
{
	REGISTER(coordinate_system, ECFMetaData()
			.setdescription({ "Material coordinate system." }));

	physical_model = PHYSICAL_MODEL::THERMAL;
	REGISTER(physical_model, ECFMetaData()
			.setdescription({ "Select used physical model." })
			.setdatatype({ ECFDataType::ENUM_FLAGS })
			.addoption(ECFOption().setname("THERMAL").setdescription("Model used by ADVECTION DIFFUSION.")
					.allowonly([&] () { return _allowed_physical_models & PHYSICAL_MODEL::THERMAL; }))
			.addoption(ECFOption().setname("LINEAR_ELASTIC").setdescription("One of models used by STRUCTURAL MECHANICS.")
					.allowonly([&] () { return _allowed_physical_models & PHYSICAL_MODEL::LINEAR_ELASTIC; })));

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

espreso::MaterialConfiguration::MaterialConfiguration(PHYSICAL_MODEL allowedPhysicalModels, bool is3D)
: espreso::MaterialConfiguration()
{
	_allowed_physical_models = allowedPhysicalModels;
	coordinate_system.is3D = is3D;
	thermal_conductivity.is3D = is3D;
	linear_elastic_properties.is3D = is3D;
}
