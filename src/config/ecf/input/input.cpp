
#include "input.h"

#include "../../configuration.hpp"

using namespace espreso;

InputConfiguration::InputConfiguration()
{
	path = ".";
	REGISTER(path, ECFMetaData()
            .setdescription({ "Path" })
			.setdatatype({ ECFDataType::STRING }));

	keep_material_sets = false;
	REGISTER(keep_material_sets, ECFMetaData()
            .setdescription({ "Keep material sets" })
			.setdatatype({ ECFDataType::BOOL }));

	convert_database = false;
	REGISTER(convert_database, ECFMetaData()
            .setdescription({ "Convert database" })
			.setdatatype({ ECFDataType::BOOL }));

	scale_factor = 1;
	REGISTER(scale_factor, ECFMetaData()
            .setdescription({ "Scale factor" })
			.setdatatype({ ECFDataType::FLOAT }));
}




