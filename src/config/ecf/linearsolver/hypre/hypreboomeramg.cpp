
#include "../../../configuration.hpp"
#include "hypreboomeramg.h"

using namespace espreso;

HYPREBoomerAMGConfiguration::HYPREBoomerAMGConfiguration()
{
	convergence_tolerance = 1e-7;;
	REGISTER(convergence_tolerance, ECFMetaData()
			.setdescription({ "Convergence tolerance" })
			.setdatatype({ ECFDataType::FLOAT }));

	min_iterations = 1;
	REGISTER(min_iterations, ECFMetaData()
			.setdescription({ "Minimum number of iterations." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	max_iterations = 1;
	REGISTER(max_iterations, ECFMetaData()
			.setdescription({ "Maximum number of iterations." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	cycle_type = CYCLE_TYPE::V_CYCLE;
	REGISTER(cycle_type, ECFMetaData()
			.setdescription({ "Cycle type...." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("V_CYCLE").setdescription("doc..."))
			.addoption(ECFOption().setname("W_CYCLE").setdescription("doc...")));
}





