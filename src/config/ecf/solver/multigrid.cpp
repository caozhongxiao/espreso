
#include "multigrid.h"
#include "config/configuration.hpp"

espreso::MultigridConfiguration::MultigridConfiguration()
{
	solver = SOLVER::CG;
	REGISTER(solver, ECFMetaData()
			.setdescription({ "HYPRE solver" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CG").setdescription("CG solver."))
			.addoption(ECFOption().setname("GMRES").setdescription("GMRES solver."))
			.addoption(ECFOption().setname("FGMRES").setdescription("FGMRES solver."))
			.addoption(ECFOption().setname("BICGS").setdescription("BICGS solver."))
			.addoption(ECFOption().setname("BICGSTAB").setdescription("BICGSTAB solver."))
			.addoption(ECFOption().setname("TFQMR").setdescription("TFQMR solver."))
			.addoption(ECFOption().setname("SYMQMR").setdescription("SYMQMR solver."))
			.addoption(ECFOption().setname("SUPERLU").setdescription("SUPERLU solver."))
			.addoption(ECFOption().setname("SUPERLUX").setdescription("SUPERLUX solver.")));

	preconditioner = PRECONDITIONER::BOOMERAMG;
	REGISTER(preconditioner, ECFMetaData()
			.setdescription({ "Preconditioner" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DIAGONAL").setdescription("DIAGONAL preconditioner."))
			.addoption(ECFOption().setname("PILUT").setdescription("PILUT preconditioner."))
			.addoption(ECFOption().setname("EUCLID").setdescription("EUCLID preconditioner."))
			.addoption(ECFOption().setname("PARASAILS").setdescription("PARASAILS preconditioner."))
			.addoption(ECFOption().setname("BOOMERAMG").setdescription("BOOMERAMG preconditioner."))
			.addoption(ECFOption().setname("POLY").setdescription("POLY preconditioner."))
			.addoption(ECFOption().setname("MLI").setdescription("MLI preconditioner.")));

	precision = 1e-5;
	REGISTER(precision, ECFMetaData()
			.setdescription({ "Precision" })
			.setdatatype({ ECFDataType::FLOAT }));

	max_iterations = 200;
	REGISTER(max_iterations, ECFMetaData()
			.setdescription({ "Max iterations" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
}



