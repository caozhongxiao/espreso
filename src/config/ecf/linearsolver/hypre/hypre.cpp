
#include "hypre.h"

#include "../../../configuration.hpp"

using namespace espreso;

HypreConfiguration::HypreConfiguration()
{
	solver_type = SOLVER_TYPE::BoomerAMG;
	REGISTER(solver_type, ECFMetaData()
			.setdescription({ "Solver type...." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("BoomerAMG").setdescription("doc..."))
			.addoption(ECFOption().setname("PCG").setdescription("doc..."))
			.addoption(ECFOption().setname("GMRES").setdescription("doc...")));

	REGISTER(boomeramg, ECFMetaData()
		.setdescription({ "BoomerAMG settings." }));

	REGISTER(pcg, ECFMetaData()
		.setdescription({ "PCG settings." }));
}



