
#include "config/configuration.hpp"
#include "heattransfer.h"

espreso::ConvectionConfiguration::ConvectionConfiguration()
: heat_transfer_coefficient(ECFMetaData::getboundaryconditionvariables()),
  external_temperature(ECFMetaData::getboundaryconditionvariables()),
  wall_height(ECFMetaData::getboundaryconditionvariables()),
  tilt_angle(ECFMetaData::getboundaryconditionvariables()),
  diameter(ECFMetaData::getboundaryconditionvariables()),
  plate_length(ECFMetaData::getboundaryconditionvariables()),
  fluid_velocity(ECFMetaData::getboundaryconditionvariables()),
  plate_distance(ECFMetaData::getboundaryconditionvariables()),
  length(ECFMetaData::getboundaryconditionvariables()),
  absolute_pressure(ECFMetaData::getboundaryconditionvariables())
{
	type = TYPE::USER;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Convection type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("USER").setdescription("User defined."))
			.addoption(ECFOption().setname("EXTERNAL_NATURAL").setdescription("External natural."))
			.addoption(ECFOption().setname("INTERNAL_NATURAL").setdescription("Internal natural."))
			.addoption(ECFOption().setname("EXTERNAL_FORCED").setdescription("External forced."))
			.addoption(ECFOption().setname("INTERNAL_FORCED").setdescription("Internal forced.")));

	variant = VARIANT::VERTICAL_WALL;
	REGISTER(variant, ECFMetaData()
			.setdescription({ "Variant" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("VERTICAL_WALL").setdescription("Vertical wall."))
			.addoption(ECFOption().setname("INCLINED_WALL").setdescription("Inclined wall."))
			.addoption(ECFOption().setname("HORIZONTAL_CYLINDER").setdescription("Horizontal cylinder."))
			.addoption(ECFOption().setname("SPHERE").setdescription("Sphere."))
			.addoption(ECFOption().setname("HORIZONTAL_PLATE_UP").setdescription("Horizontal place up."))
			.addoption(ECFOption().setname("HORIZONTAL_PLATE_DOWN").setdescription("Horizontal plate down."))
			.addoption(ECFOption().setname("AVERAGE_PLATE").setdescription("Average plate."))
			.addoption(ECFOption().setname("PARALLEL_PLATES").setdescription("Parallel plates."))
			.addoption(ECFOption().setname("CIRCULAR_TUBE").setdescription("Circular tube."))
			.addoption(ECFOption().setname("TUBE").setdescription("Tube.")));

	fluid = FLUID::AIR;
	REGISTER(fluid, ECFMetaData()
			.setdescription({ "Fluid type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("AIR").setdescription("Air."))
			.addoption(ECFOption().setname("WATER").setdescription("Water."))
			.addoption(ECFOption().setname("ENGINE_OIL").setdescription("Engine oil."))
			.addoption(ECFOption().setname("TRANSFORMER_OIL").setdescription("Tranformer oil.")));

	REGISTER(heat_transfer_coefficient, ECFMetaData()
			.setdescription({ "Heat transfer coefficient" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(external_temperature, ECFMetaData()
			.setdescription({ "Ambient temperature" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(wall_height, ECFMetaData()
			.setdescription({ "Wall height" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(tilt_angle, ECFMetaData()
			.setdescription({ "Tilt angle" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(diameter, ECFMetaData()
			.setdescription({ "Diameter" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(plate_length, ECFMetaData()
			.setdescription({ "Plate length" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(fluid_velocity, ECFMetaData()
			.setdescription({ "Fluid velocity" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(plate_distance, ECFMetaData()
			.setdescription({ "Plate distance" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(length, ECFMetaData()
			.setdescription({ "Length" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(absolute_pressure, ECFMetaData()
			.setdescription({ "Absolute pressure" })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

espreso::RadiationConfiguration::RadiationConfiguration()
: emissivity(ECFMetaData::getboundaryconditionvariables()),
  external_temperature(ECFMetaData::getboundaryconditionvariables())
{
	REGISTER(emissivity, ECFMetaData()
			.setdescription({ "Emissivity" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(external_temperature, ECFMetaData()
			.setdescription({ "Ambient temperature" })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

espreso::HeatTransferLoadStepConfiguration::HeatTransferLoadStepConfiguration(DIMENSION dimension)
: LoadStepConfiguration("temperature", "heat")
{
	REGISTER(temperature, ECFMetaData()
			.setdescription({ "The name of a region.", "Temperature" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "273.15" }),
			ECFMetaData::getboundaryconditionvariables());
	REGISTER(heat_source, ECFMetaData()
			.setdescription({ "The name of a region.", "Heat source" })
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "273.15" }),
			ECFMetaData::getboundaryconditionvariables());
	REGISTER(translation_motions, ECFMetaData()
			.setdescription({ "The name of a region.", "Translation motion" })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION" }),
			dimension, ECFMetaData::getboundaryconditionvariables(), "0");

	REGISTER(heat_flux, ECFMetaData()
			.setdescription({ "The name of a region.", "Heat flux" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "500" }),
			ECFMetaData::getboundaryconditionvariables());
	REGISTER(heat_flow, ECFMetaData()
			.setdescription({ "The name of a region.", "Heat flow" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "500" }),
			ECFMetaData::getboundaryconditionvariables());

	REGISTER(convection, ECFMetaData()
			.setdescription({ "The name of a region.", "Convection" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" }));
	REGISTER(diffuse_radiation, ECFMetaData()
			.setdescription({ "The name of a region.", "Diffuse radiation" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" }));
}

espreso::HeatTransferConfiguration::HeatTransferConfiguration(DIMENSION dimension)
: PhysicsConfiguration(dimension, MaterialConfiguration::PHYSICAL_MODEL::THERMAL)
{
	REGISTER(materials, ECFMetaData()
			.setdescription({ "The name of a material.", "Material description" })
			.setdatatype({ ECFDataType::STRING })
			.setpattern({ "MY_MATERIAL" }),
			dimension, MaterialConfiguration::PHYSICAL_MODEL::THERMAL);
	ecfdescription->moveLastBefore(PNAME(material_set));

	stabilization = STABILIZATION::SUPG;
	REGISTER(stabilization, ECFMetaData()
			.setdescription({ "Inconsistent stabilization" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("SUPG").setdescription("SUPG stabilization"))
			.addoption(ECFOption().setname("CAU").setdescription("CAU stabilization")));

	sigma = 0;
	REGISTER(sigma, ECFMetaData()
			.setdescription({ "Inconsistent stabilization parameter" })
			.setdatatype({ ECFDataType::FLOAT }));

	init_temp_respect_bc = true;
	REGISTER(init_temp_respect_bc, ECFMetaData()
			.setdescription({ "Initial temperature follows BC" })
			.setdatatype({ ECFDataType::BOOL }));

	diffusion_split = false;
	REGISTER(diffusion_split, ECFMetaData()
			.setdescription({ "Thermal shock stabilization" })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(initial_temperature, ECFMetaData()
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
			.setdescription({ "The name of a region.", "Initial temperature" })
			.setpattern({ "MY_REGION", "273.15" }),
			ECFMetaData().getboundaryconditionvariables(), "273.15");

	if (dimension == DIMENSION::D2) {
		REGISTER(thickness, ECFMetaData()
				.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
				.setdescription({ "The name of a region.", "Thickness" })
				.setpattern({ "MY_REGION", "1" }),
				ECFMetaData().getboundaryconditionvariables(), "1");
	}

	REGISTER(load_steps_settings, ECFMetaData()
			.setdescription({ "Settings for each load step", "LoadStep" })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" }),
			dimension);
}


