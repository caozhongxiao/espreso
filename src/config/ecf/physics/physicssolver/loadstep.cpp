
#include "loadstep.h"
#include "config/configuration.hpp"

espreso::LoadStepSolverConfiguration::LoadStepSolverConfiguration()
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
			.addoption(ECFOption().setname("HYPRE").setdescription("Use hypre library."))
			.addoption(ECFOption().setname("MKLPDSS").setdescription("Use parallel direct sparse solver from MKL.")));


	REGISTER(feti, ECFMetaData()
			.setdescription({ "FETI solver settings" })
			.allowonly([&] () { return solver == SOLVER::FETI; }));
	REGISTER(hypre, ECFMetaData()
			.setdescription({ "HYPRE multigrid solver settings" })
			.allowonly([&] () { return solver == SOLVER::HYPRE; }));
	REGISTER(mklpdss, ECFMetaData()
			.setdescription({ "MKL parallel direct sparse solver" })
			.allowonly([&] () { return solver == SOLVER::MKLPDSS; }));
}

espreso::HeatTransferLoadStepSolverConfiguration::HeatTransferLoadStepSolverConfiguration()
: nonlinear_solver("temperature", "heat")
{
	REGISTER(nonlinear_solver, ECFMetaData()
			.setdescription({ "Non-linear physics solver settings" })
			.allowonly([&] () { return mode == MODE::NONLINEAR; }));

	REGISTER(transient_solver, ECFMetaData()
			.setdescription({ "Transient physics solver settings" })
			.allowonly([&] () { return type == TYPE::TRANSIENT; }));
}

espreso::StructuralMechanicsLoadStepSolverConfiguration::StructuralMechanicsLoadStepSolverConfiguration()
: nonlinear_solver("displacement", "forces")
{
	REGISTER(nonlinear_solver, ECFMetaData()
			.setdescription({ "Non-linear physics solver settings" })
			.allowonly([&] () { return mode == MODE::NONLINEAR; }));

	REGISTER(transient_solver, ECFMetaData()
			.setdescription({ "Transient physics solver settings" })
			.allowonly([&] () { return type == TYPE::TRANSIENT; }));
}


