
#include "transientsolver.h"
#include "../../configuration.hpp"

espreso::AutoTimeSteppingConfiguration::AutoTimeSteppingConfiguration()
{
	allowed = true;
	REGISTER(allowed, ECFMetaData()
			.setdescription({ "Allow auto time stepping." })
			.setdatatype({ ECFDataType::BOOL }));

	initial_time_step = 1e-2;
	min_time_step = 1e-3;
	max_time_step = 1;
	REGISTER(initial_time_step, ECFMetaData()
			.setdescription({ "Initial time step." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(min_time_step, ECFMetaData()
			.setdescription({ "Minimal time step." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(max_time_step, ECFMetaData()
			.setdescription({ "Maximal time step." })
			.setdatatype({ ECFDataType::FLOAT }));
}

espreso::TransientSolverConfiguration::TransientSolverConfiguration()
{
	method = METHOD::CRANK_NICOLSON;
	REGISTER(method, ECFMetaData()
			.setdescription({ "Transien solver used method." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CRANK_NICOLSON").setdescription("Alpha = 0.5."))
			.addoption(ECFOption().setname("FORWARD_DIFF").setdescription("Alpha = ??."))
			.addoption(ECFOption().setname("GALERKIN").setdescription("Alpha = 2 / 3."))
			.addoption(ECFOption().setname("BACKWARD_DIFF").setdescription("Alpha = 1."))
			.addoption(ECFOption().setname("USER").setdescription("User defined Alpha from interval <0, 1).")));

	alpha = 0.5;
	time_step = 0.1;
	REGISTER(alpha, ECFMetaData()
			.setdescription({ "User defined Alpha." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(time_step, ECFMetaData()
			.setdescription({ "Duration of a time step." })
			.setdatatype({ ECFDataType::FLOAT }));
}



