
#include "transientsolver.h"
#include "../../../configuration.hpp"

espreso::AutoTimeSteppingConfiguration::AutoTimeSteppingConfiguration()
{
	allowed = false;
	REGISTER(allowed, ECFMetaData()
			.setdescription({ "Allow auto time stepping." })
			.setdatatype({ ECFDataType::BOOL }));

	addSeparator();

	min_time_step = 1e-3;
	max_time_step = 1;
	REGISTER(min_time_step, ECFMetaData()
			.setdescription({ "Minimal time step." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(max_time_step, ECFMetaData()
			.setdescription({ "Maximal time step." })
			.setdatatype({ ECFDataType::FLOAT }));

	addSeparator();

	oscilation_limit = 0.5;
	IDFactor = 3;
	REGISTER(oscilation_limit, ECFMetaData()
			.setdescription({ "Oscilation limit." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(IDFactor, ECFMetaData()
			.setdescription({ "Increace/ decrease factor." })
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
	REGISTER(alpha, ECFMetaData()
			.setdescription({ "User defined Alpha." })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return method == METHOD::USER; }));

	addSpace();

	time_step = 0.1;
	REGISTER(time_step, ECFMetaData()
				.setdescription({ "Duration of a time step." })
				.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(auto_time_stepping, ECFMetaData()
			.setdescription({ "Auto time stepping." }));


}



