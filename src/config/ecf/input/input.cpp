
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

	convert_database = false;
	REGISTER(convert_database, ECFMetaData()
			.setdescription({ "Read data, create visual output and store to ESPRESO binary format." })
			.setdatatype({ ECFDataType::BOOL }));
}




