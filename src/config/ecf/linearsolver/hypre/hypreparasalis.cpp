
#include "hypreparasalis.h"

#include "config/configuration.hpp"

using namespace espreso;

HYPREParasalisConfiguration::HYPREParasalisConfiguration()
{
	convergence_tolerance = 1e-8;
	REGISTER(convergence_tolerance, ECFMetaData()
			.setdescription({ "Set the convergence tolerance" })
			.setdatatype({ ECFDataType::FLOAT }));

}
