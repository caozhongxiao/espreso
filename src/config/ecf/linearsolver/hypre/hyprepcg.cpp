
#include "hyprepcg.h"

#include "../../../configuration.hpp"

using namespace espreso;

HYPREPCGConfiguration::HYPREPCGConfiguration()
{
	relative_conv_tol = 1e-8;
	REGISTER(relative_conv_tol, ECFMetaData()
			.setdescription({ "Relative convergence tolerance" })
			.setdatatype({ ECFDataType::FLOAT }));

	preconditioner = PRECONDITIONER::BoomerAMG;
	REGISTER(preconditioner, ECFMetaData()
			.setdescription({ "Preconditioner" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("BoomerAMG").setdescription("doc..."))
			.addoption(ECFOption().setname("Parasalis").setdescription("doc...")));

	REGISTER(boomeramg, ECFMetaData()
			.setdescription({ "BoomerAMG settings." }));
}


