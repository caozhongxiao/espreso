
#include "solver.h"
#include "../../configuration.hpp"

espreso::PhysicsSolverConfiguration::PhysicsSolverConfiguration()
{
	interpolation = INTERPOLATION::LINEAR;
	REGISTER(interpolation, ECFMetaData()
			.setdescription({ "Assembler data interpolation type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR").setdescription("Linear interpolation."))
			.addoption(ECFOption().setname("QUADRATIC").setdescription("Quadratic interpolation.")));

	load_steps = 1;
	REGISTER(load_steps, ECFMetaData()
			.setdescription({ "Number of load steps." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
}

espreso::AdvectionDiffusionPhysicsSolverConfiguration::AdvectionDiffusionPhysicsSolverConfiguration()
{
	REGISTER(load_steps_settings, ECFMetaData()
			.setdescription({ "Load step index", "Load steps settings." })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" }));
}

espreso::StructuralMechanicsPhysicsSolverConfiguration::StructuralMechanicsPhysicsSolverConfiguration()
{
	REGISTER(load_steps_settings, ECFMetaData()
			.setdescription({ "Load step index", "Load steps settings." })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" }));
}

