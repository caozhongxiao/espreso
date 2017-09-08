
#include "structuralmechanics.h"
#include "../../configuration.hpp"

espreso::StructuralMechanicsLoadStepConfiguration::StructuralMechanicsLoadStepConfiguration(DIMENSION dimension)
: LoadStepConfiguration("displacement", "forces")
{
	REGISTER(temperature, ECFMetaData()
			.setdescription({ "The name of a region.", "Temperature of a given region." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "273.15" }));
	REGISTER(displacement, ECFMetaData()
			.setdescription({ "The name of a region.", "Fixed displacement of a given region." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", dimension == DIMENSION::D2 ? "X 0, Y 0" : "X 0, Y 0, Z 0" }));
	REGISTER(acceleration, ECFMetaData()
			.setdescription({ "The name of a region.", "Acceleration of a given region." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", dimension == DIMENSION::D2 ? "X 0, Y -9.81" : "X 0, Y 0, Z -9.81" }));
	REGISTER(normal_pressure, ECFMetaData()
			.setdescription({ "The name of a region.", "Normal pressure on a given region." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "0" }));
	REGISTER(obstacle, ECFMetaData()
			.setdescription({ "The name of a region.", "Obstacle for a given region." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "2" }));
	REGISTER(normal_direction, ECFMetaData()
			.setdescription({ "The name of a region.", "Normal direction of a given region." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "x" }));
}

espreso::StructuralMechanicsConfiguration::StructuralMechanicsConfiguration(DIMENSION dimension)
: PhysicsConfiguration(dimension), element_behaviour(ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS)
{
	REGISTER(
			materials,
			ECFMetaData()
				.setdescription({ "The name of a material.", "Material description." })
				.setdatatype({ ECFDataType::STRING })
				.setpattern({ "MY_MATERIAL" }),
			dimension, MaterialConfiguration::PHYSICAL_MODEL::LINEAR_ELASTIC);
	moveLastBefore(PNAME(material_set));

	REGISTER(element_behaviour, ECFMetaData()
			.setdescription({ "Physics solver type." })
			.setdatatype({ ECFDataType::OPTION })
			.allowonly([&] () { return dimension == DIMENSION::D2; })
			.addoption(ECFOption().setname("PLANE_STRAIN").setdescription("Plane strain."))
			.addoption(ECFOption().setname("AXISYMMETRIC").setdescription("Axisymmetric."))
			.addoption(ECFOption().setname("PLANE_STRESS").setdescription("Plane stress."))
			.addoption(ECFOption().setname("PLANE_STRESS_WITH_THICKNESS").setdescription("Plane stress with thickness.")));
	moveLastBefore(PNAME(materials));

	REGISTER(
			load_steps_settings,
			ECFMetaData()
				.setdescription({ "Settings for each load step.", "Load step index." })
				.setdatatype({ ECFDataType::LOAD_STEP })
				.setpattern({ "1" }),
			dimension);
}




