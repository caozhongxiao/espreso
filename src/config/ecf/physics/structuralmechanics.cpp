
#include "structuralmechanics.h"
#include "../../configuration.hpp"

espreso::StructuralMechanicsConfiguration::StructuralMechanicsConfiguration()
{
	REGISTER(physics_solver, ECFMetaData()
			.setdescription({ "Physics solver settings." }));

	REGISTER(initial_temperature, ECFMetaData()
			.setdescription({ "The name of a region.", "Initial temperature of a given region." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "273.15" }));

	REGISTER(temperature, ECFMetaData()
			.setdescription({ "Load step index.", "The name of a region.", "Temperature of a given region." })
			.setdatatype({ ECFDataType::LOAD_STEP, ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "1", "MY_REGION", "273.15" }));
	REGISTER(displacement, ECFMetaData()
			.setdescription({ "Load step index.", "The name of a region.", "Fixed displacement of a given region." })
			.setdatatype({ ECFDataType::LOAD_STEP, ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "1", "MY_REGION", "X 0, Y 0, Z 0" }));
	REGISTER(acceleration, ECFMetaData()
			.setdescription({ "Load step index.", "The name of a region.", "Acceleration of a given region." })
			.setdatatype({ ECFDataType::LOAD_STEP, ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "1", "MY_REGION", "X 0, Y 0, Z -9.81" }));
	REGISTER(normal_pressure, ECFMetaData()
			.setdescription({ "Load step index.", "The name of a region.", "Normal pressure on a given region." })
			.setdatatype({ ECFDataType::LOAD_STEP, ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "1", "MY_REGION", "0" }));
	REGISTER(obstacle, ECFMetaData()
			.setdescription({ "Load step index.", "The name of a region.", "Obstacle for a given region." })
			.setdatatype({ ECFDataType::LOAD_STEP, ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "1", "MY_REGION", "2" }));
	REGISTER(normal_direction, ECFMetaData()
			.setdescription({ "Load step index.", "The name of a region.", "Normal direction of a given region." })
			.setdatatype({ ECFDataType::LOAD_STEP, ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "1", "MY_REGION", "x" }));

	REGISTER(material_set, ECFMetaData()
			.setdescription({ "The name of a region.", "The name of a material." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::MATERIAL })
			.setpattern({ "MY_REGION", "MY_MATERIAL" }));

	post_process = false;
	REGISTER(post_process, ECFMetaData()
			.setdescription({ "Turn on/off post-process." })
			.setdatatype({ ECFDataType::BOOL}));
}

espreso::StructuralMechanics2DMaterialConfiguration::StructuralMechanics2DMaterialConfiguration()
{
	physical_model = PHYSICAL_MODEL::LINEAR_ELASTIC;
	coordinationSystem.is3D = false;
	linear_elastic_properties.is3D = false;
}

espreso::StructuralMechanics3DMaterialConfiguration::StructuralMechanics3DMaterialConfiguration()
{
	physical_model = PHYSICAL_MODEL::LINEAR_ELASTIC;
}

espreso::StructuralMechanics2DConfiguration::StructuralMechanics2DConfiguration()
{
	getWithError(PNAME(displacement))->metadata.pattern[2] = "X 0, Y 0";
	getWithError(PNAME(acceleration))->metadata.pattern[2] = "X 0, Y -9.81";
	REGISTER(thickness, ECFMetaData()
			.setdescription({ "Load step index.", "The name of a region.", "Thickness of a given region." })
			.setdatatype({ ECFDataType::LOAD_STEP, ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setpattern({ "1", "MY_REGION", "1" }));
	moveLastBefore(PNAME(initial_temperature));

	REGISTER(materials, ECFMetaData()
			.setdescription({ "The name of a material.", "Material description." })
			.setdatatype({ ECFDataType::STRING })
			.setpattern({ "MY_MATERIAL" }));
	moveLastBefore(PNAME(material_set));

	element_behaviour = ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS;
	REGISTER(element_behaviour, ECFMetaData()
			.setdescription({ "Physics solver type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("PLANE_STRAIN").setdescription("Plane strain."))
			.addoption(ECFOption().setname("AXISYMMETRIC").setdescription("Axisymmetric."))
			.addoption(ECFOption().setname("PLANE_STRESS").setdescription("Plane stress."))
			.addoption(ECFOption().setname("PLANE_STRESS_WITH_THICKNESS").setdescription("Plane stress with thickness.")));
	moveLastBefore(PNAME(thickness));
}

espreso::StructuralMechanics3DConfiguration::StructuralMechanics3DConfiguration()
{
	REGISTER(materials, ECFMetaData()
			.setdescription({ "The name of a material.", "Material description." })
			.setdatatype({ ECFDataType::STRING })
			.setpattern({ "MY_MATERIAL" }));
	moveLastBefore(PNAME(material_set));

	discretization = DISCRETIZATION::FEM;
	REGISTER(discretization, ECFMetaData()
			.setdescription({ "Discretization of stiffness matrices." })
		.setdatatype({ ECFDataType::OPTION })
		.addoption(ECFOption().setname("FEM").setdescription("Finite elements."))
		.addoption(ECFOption().setname("BEM").setdescription("Boundary elements.")));
}




