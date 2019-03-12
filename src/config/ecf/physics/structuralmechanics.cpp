
#include "structuralmechanics.h"
#include "config/configuration.hpp"

espreso::StructuralMechanicsLoadStepConfiguration::StructuralMechanicsLoadStepConfiguration(DIMENSION *dimension)
{
	large_displacement = false;
	REGISTER(large_displacement, ECFMetaData()
			.setdescription({ "Turn on large displacement." })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(temperature, ECFMetaData()
			.setdescription({ "The name of a region.", "Temperature of a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "275.15" }),
			ECFMetaData::getboundaryconditionvariables());
	REGISTER(normal_pressure, ECFMetaData()
			.setdescription({ "The name of a region.", "Normal pressure on a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "0" }),
			ECFMetaData::getboundaryconditionvariables());

	REGISTER(angular_velocity, ECFMetaData()
			.setdescription({ "The name of a region.", "Angular velocity of a given region." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION" }),
			dimension, ECFMetaData::getboundaryconditionvariables());

	REGISTER(normal_direction, ECFMetaData()
			.setdescription({ "The name of a region.", "Normal direction of a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" }),
			dimension, ECFMetaData::getboundaryconditionvariables());

	REGISTER(obstacle, ECFMetaData()
			.setdescription({ "The name of a region.", "Obstacle for a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" }),
			dimension, ECFMetaData::getboundaryconditionvariables());

	REGISTER(acceleration, ECFMetaData()
			.setdescription({ "The name of a region.", "Acceleration of a given region." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION" }),
			dimension, ECFMetaData::getboundaryconditionvariables(), "0");


	REGISTER(displacement, ECFMetaData()
			.setdescription({ "The name of a region.", "Fixed displacement of a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" }),
			dimension, ECFMetaData::getboundaryconditionvariables());
}

espreso::StructuralMechanicsConfiguration::StructuralMechanicsConfiguration(DIMENSION D)
: PhysicsConfiguration(D, MaterialConfiguration::PHYSICAL_MODEL::STRUCTURAL_MECHANICS)
{
	REGISTER(materials, ECFMetaData()
				.setdescription({ "The name of a material.", "Material description." })
				.setdatatype({ ECFDataType::STRING })
				.setpattern({ "MY_MATERIAL" }),
				&dimension, MaterialConfiguration::PHYSICAL_MODEL::STRUCTURAL_MECHANICS);
	ecfdescription->moveLastBefore(PNAME(material_set));

	element_behaviour = ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS;
	REGISTER(element_behaviour, ECFMetaData()
			.setdescription({ "Physics solver type." })
			.setdatatype({ ECFDataType::OPTION })
			.allowonly([&] () { return this->dimension == DIMENSION::D2; })
			.addoption(ECFOption().setname("PLANE_STRAIN").setdescription("Plane strain."))
			.addoption(ECFOption().setname("AXISYMMETRIC").setdescription("Axisymmetric."))
			.addoption(ECFOption().setname("PLANE_STRESS").setdescription("Plane stress."))
			.addoption(ECFOption().setname("PLANE_STRESS_WITH_THICKNESS").setdescription("Plane stress with thickness.")));
	ecfdescription->moveLastBefore(PNAME(materials));

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
			.setdescription({ "Settings for each load step.", "Load step index." })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" }),
			&dimension);
}




