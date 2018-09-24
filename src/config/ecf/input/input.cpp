
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

	addSeparator();

	REGISTER(load_by, ECFMetaData()
			.setdescription({ "Load by. " })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NODES").setdescription("Computation nodes."))
			.addoption(ECFOption().setname("PROCESSES").setdescription("MPI processes.")));

	REGISTER(load_pattern, ECFMetaData()
			.setdescription({ "Load pattern. " })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("PREFIX").setdescription("Only the first subset loads data."))
			.addoption(ECFOption().setname("SUBSET").setdescription("Each n-th node(process) loads data.")));

	REGISTER(pattern_value, ECFMetaData()
			.setdescription({ "Define the fist n-process / each n-th process that loads data." })
			.setdatatype({ ECFDataType::INTEGER }));
}




