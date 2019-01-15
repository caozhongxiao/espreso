
#include "input.h"

#include "config/configuration.hpp"

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

	ecfdescription->addSeparator();

	granularity = ProcessesReduction::Granularity::PROCESSES;
	REGISTER(granularity, ECFMetaData()
			.setdescription({ "A granularity of loaders. " })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NODES").setdescription("Computation nodes."))
			.addoption(ECFOption().setname("PROCESSES").setdescription("MPI processes.")));

	pattern = ProcessesReduction::Pattern::SUBSET;
//	REGISTER(pattern, ECFMetaData()
//			.setdescription({ "A pattern for grouping processes. " })
//			.setdatatype({ ECFDataType::OPTION })
//			.addoption(ECFOption().setname("PREFIX").setdescription("Only the first subset loads data."))
//			.addoption(ECFOption().setname("SUBSET").setdescription("Each n-th node(process) loads data.")));

	reduction_ratio = 1;
	REGISTER(reduction_ratio, ECFMetaData()
			.setdescription({ "Defines the reduction ratio." })
			.setdatatype({ ECFDataType::INTEGER }));
}




