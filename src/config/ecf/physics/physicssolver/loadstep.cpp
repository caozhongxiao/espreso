
#include "loadstep.h"
#include "../../../configuration.hpp"

espreso::LoadStepConfiguration::LoadStepConfiguration(const std::string &firstResidualName, const std::string &secondResidualName)
: nonlinear_solver(firstResidualName, secondResidualName)
{
	duration_time = 1;
	REGISTER(duration_time, ECFMetaData()
			.setdescription({ "Duration of a load step." })
			.setdatatype({ ECFDataType::FLOAT }));

	type = TYPE::STEADY_STATE;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Physics solver type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("STEADY_STATE").setdescription("Steady state load step."))
			.addoption(ECFOption().setname("TRANSIENT").setdescription("Transient load step.")));

	mode = MODE::LINEAR;
	REGISTER(mode, ECFMetaData()
			.setdescription({ "Physics solver mode - set according to material properties." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR").setdescription("Linear material behavior."))
			.addoption(ECFOption().setname("NONLINEAR").setdescription("Nonlinear material behavior.")));

	solver = SOLVER::FETI;
	REGISTER(solver, ECFMetaData()
			.setdescription({ "Used linear solver method." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("FETI").setdescription("Use ESPRESO as linear solver."))
			.addoption(ECFOption().setname("MULTIGRID").setdescription("Use hypre library as MULTIGRID solver.")));

	REGISTER(nonlinear_solver, ECFMetaData()
				.setdescription({ "Non-linear physics solver settings." })
				.allowonly([&] () { return mode == MODE::NONLINEAR; }));
	REGISTER(transient_solver, ECFMetaData()
			.setdescription({ "Transient physics solver settings." })
			.allowonly([&] () { return type == TYPE::TRANSIENT; }));

	REGISTER(feti, ECFMetaData()
			.setdescription({ "ESPRESO FETI solver settings." }));
	REGISTER(multigrid, ECFMetaData()
			.setdescription({ "Hypre multigrid solver settings." }));
}


