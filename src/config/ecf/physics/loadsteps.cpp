
#include "loadsteps.h"
#include "../../configuration.hpp"

espreso::LoadStepsConfiguration::LoadStepsConfiguration()
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

	REGISTER(transient_solver, ECFMetaData()
			.setdescription({ "Transient physics solver settings." }));


	REGISTER(feti, ECFMetaData()
			.setdescription({ "ESPRESO FETI solver settings." }));
	REGISTER(multigrid, ECFMetaData()
			.setdescription({ "Hypre multigrid solver settings." }));
}

espreso::AdvectionDiffusionLoadStepsConfiguration::AdvectionDiffusionLoadStepsConfiguration()
{
	REGISTER(nonlinear_solver, ECFMetaData()
			.setdescription({ "Non-linear physics solver settings." }));
	moveLastBefore(PNAME(transient_solver));
}

espreso::StructuralMechanicsLoadStepsConfiguration::StructuralMechanicsLoadStepsConfiguration()
{
	REGISTER(nonlinear_solver, ECFMetaData()
			.setdescription({ "Non-linear physics solver settings." }));
	moveLastBefore(PNAME(transient_solver));
}



