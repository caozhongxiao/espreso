
#include "loadstep.h"
#include "config/configuration.hpp"

espreso::LoadStepConfiguration::LoadStepConfiguration(const std::string &firstResidualName, const std::string &secondResidualName)
: nonlinear_solver(firstResidualName, secondResidualName)
{
	duration_time = 1;
	REGISTER(duration_time, ECFMetaData()
			.setdescription({ "Duration" })
			.setdatatype({ ECFDataType::FLOAT }));

	type = TYPE::STEADY_STATE;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Simulation type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("STEADY_STATE").setdescription("Steady state load step."))
			.addoption(ECFOption().setname("TRANSIENT").setdescription("Transient load step.")));

	mode = MODE::LINEAR;
	REGISTER(mode, ECFMetaData()
			.setdescription({ "Simulation mode" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR").setdescription("Linear material behavior."))
			.addoption(ECFOption().setname("NONLINEAR").setdescription("Nonlinear material behavior.")));

	solver = SOLVER::FETI;
	REGISTER(solver, ECFMetaData()
			.setdescription({ "Linear solver" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("FETI").setdescription("Use ESPRESO as linear solver."))
			.addoption(ECFOption().setname("MULTIGRID").setdescription("Use hypre library as MULTIGRID solver."))
			.addoption(ECFOption().setname("HYPRE").setdescription("Use hypre library.")));

	REGISTER(nonlinear_solver, ECFMetaData()
			.setdescription({ "Non-linear physics solver settings" })
			.allowonly([&] () { return mode == MODE::NONLINEAR; }));
	REGISTER(transient_solver, ECFMetaData()
			.setdescription({ "Transient physics solver settings" })
			.allowonly([&] () { return type == TYPE::TRANSIENT; }));

	REGISTER(feti, ECFMetaData()
			.setdescription({ "FETI solver settings" })
			.allowonly([&] () { return solver == SOLVER::FETI; }));
	REGISTER(multigrid, ECFMetaData()
			.setdescription({ "HYPRE multigrid solver settings" })
			.allowonly([&] () { return solver == SOLVER::MULTIGRID; }));
	REGISTER(hypre, ECFMetaData()
			.setdescription({ "HYPRE multigrid solver settings" })
			.allowonly([&] () { return solver == SOLVER::HYPRE; }));
}


