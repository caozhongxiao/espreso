
#include "physics.h"
#include "../../configuration.hpp"

espreso::PhysicsConfiguration::PhysicsConfiguration(DIMENSION dimension, MaterialConfiguration::PHYSICAL_MODEL physicalModel)
: dimension(dimension), physical_model(physicalModel)
{
	load_steps = 1;
	REGISTER(load_steps, ECFMetaData()
			.setdescription({ "Number of load steps." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	addSpace();

	interpolation = INTERPOLATION::LINEAR;
	REGISTER(interpolation, ECFMetaData()
			.setdescription({ "Data interpolation." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR").setdescription("Linear interpolation."))
			.addoption(ECFOption().setname("QUADRATIC").setdescription("Quadratic interpolation.")));

	assembler = ASSEMBLER::ELEMENTS;
	REGISTER(assembler, ECFMetaData()
			.setdescription({ "The type of assembler." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ELEMENTS").setdescription("Element based."))
			.addoption(ECFOption().setname("FACES").setdescription("Face based.")));

	REGISTER(discretization, ECFMetaData()
			.setdescription({ "Discretization settings for regions.", "Discretization of stiffness matrices." })
		.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::OPTION })
		.setpattern({ "MY_REGION", "FEM" })
		.addoption(ECFOption().setname("FEM").setdescription("Finite elements."))
		.addoption(ECFOption().setname("BEM").setdescription("Boundary elements.")));

	addSeparator();

	REGISTER(material_set, ECFMetaData()
			.setdescription({ "The name of a region.", "The name of a material." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::MATERIAL })
			.setpattern({ "MY_REGION", "MY_MATERIAL" }));

	addSeparator();

	REGISTER(initial_temperature, ECFMetaData()
			.setdescription({ "The name of a region.", "Initial temperature of a given region." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "273.15" }),
			ECFMetaData().getboundaryconditionvariables(), "273.15");

	if (dimension == DIMENSION::D2) {
		REGISTER(thickness, ECFMetaData()
				.setdescription({ "The name of a region.", "Thickness of a given region." })
				.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
				.setpattern({ "MY_REGION", "1" }),
				ECFMetaData().getboundaryconditionvariables(), "1");
	}

	addSeparator();

	contact_interfaces = false;
	REGISTER(contact_interfaces, ECFMetaData()
			.setdescription({ "Check contacts." })
			.setdatatype({ ECFDataType::BOOL }));
}



