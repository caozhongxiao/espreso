
#include "input.h"

#include "../../configuration.hpp"

using namespace espreso;

InputConfiguration::InputConfiguration()
{
	path = ".";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Path to an input example." })
			.setdatatype({ ECFDataType::STRING }));

	keep_material_sets = false;
	REGISTER(keep_material_sets, ECFMetaData()
			.setdescription({ "Keep material sets from input format." })
			.setdatatype({ ECFDataType::BOOL }));
}




