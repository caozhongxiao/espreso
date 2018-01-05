
#include "physics.h"
#include "../../configuration.hpp"

espreso::MortarConfiguration::MortarConfiguration()
{
	REGISTER(master, ECFMetaData()
			.setdescription({ "Master region." })
			.setdatatype({ ECFDataType::REGION }));

	REGISTER(slave, ECFMetaData()
			.setdescription({ "Slave region." })
			.setdatatype({ ECFDataType::REGION }));
}

espreso::ECFExpressionVector::ECFExpressionVector(DIMENSION dimension, bool fillWithZeros)
{
	REGISTER(x, ECFMetaData()
			.setdescription({ "x-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setboundaryconditionvariables());
	REGISTER(y, ECFMetaData()
			.setdescription({ "y-direction." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setboundaryconditionvariables());
	if (dimension == DIMENSION::D3) {
		REGISTER(z, ECFMetaData()
				.setdescription({ "z-direction." })
				.setdatatype({ ECFDataType::EXPRESSION })
				.setboundaryconditionvariables());
	}

	if (fillWithZeros) {
		x.value = "0";
		y.value = "0";
		x.createEvaluator(ECFMetaData().setboundaryconditionvariables().variables);
		y.createEvaluator(ECFMetaData().setboundaryconditionvariables().variables);
		if (dimension == DIMENSION::D3) {
			z.value = "0";
			z.createEvaluator(ECFMetaData().setboundaryconditionvariables().variables);
		}
	}
}

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

	REGISTER(discretization, ECFMetaData()
			.setdescription({ "Discretization settings for regions.", "Discretization of stiffness matrices." })
		.setdatatype({ ECFDataType::REGION, ECFDataType::OPTION })
		.setpattern({ "MY_REGION", "FEM" })
		.addoption(ECFOption().setname("FEM").setdescription("Finite elements."))
		.addoption(ECFOption().setname("BEM").setdescription("Boundary elements.")));

	addSeparator();

	REGISTER(material_set, ECFMetaData()
			.setdescription({ "The name of a region.", "The name of a material." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::MATERIAL })
			.setpattern({ "MY_REGION", "MY_MATERIAL" }));

	addSeparator();

	REGISTER(initial_temperature, ECFMetaData()
			.setdescription({ "The name of a region.", "Initial temperature of a given region." })
			.setdatatype({ ECFDataType::REGION, ECFDataType::EXPRESSION })
			.setcoordinatevariables()
			.setpattern({ "MY_REGION", "273.15" }));

	if (dimension == DIMENSION::D2) {
		REGISTER(thickness, ECFMetaData()
				.setdescription({ "The name of a region.", "Thickness of a given region." })
				.setdatatype({ ECFDataType::REGION, ECFDataType::EXPRESSION })
				.setcoordinatevariables()
				.setpattern({ "MY_REGION", "1" }));
	}

	addSeparator();

	REGISTER(mortar, ECFMetaData()
			.setdescription({ "Mortar interface." }));
}



