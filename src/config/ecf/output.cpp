
#include "output.h"
#include "../configuration.hpp"

espreso::OutputConfiguration::OutputConfiguration()
{
	path = "results";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Output path." })
			.setdatatype({ ECFDataType::STRING }));

	format = FORMAT::VTK_XML_ASCII;
	REGISTER(format, ECFMetaData()
			.setdescription({ "Output data format." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("VTK_LEGACY").setdescription("VTK legacy format."))
			.addoption(ECFOption().setname("VTK_XML_ASCII").setdescription("VTK XML ASCII format."))
			.addoption(ECFOption().setname("VTK_XML_BINARY").setdescription("VTK XML binary format."))
			.addoption(ECFOption().setname("ENSIGHT").setdescription("EnSight format.")));

	mode = MODE::THREAD;
	REGISTER(mode, ECFMetaData()
			.setdescription({ "ASYNC library mode." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("SYNC").setdescription("Output is synchronized."))
			.addoption(ECFOption().setname("THREAD").setdescription("Output is done by asynchronous thread."))
			.addoption(ECFOption().setname("MPI").setdescription("Output is done by separated MPI processes.")));

	output_node_group_size = 7;
	REGISTER(output_node_group_size, ECFMetaData()
			.setdescription({ "Number of MPI processes that send output data to the same storing node." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	solution = true;
	subsolution = settings = FETI_data = catalyst = false;
	catalyst_sleep_time = 0;
	REGISTER(solution, ECFMetaData()
			.setdescription({ "Store solution." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(subsolution, ECFMetaData()
			.setdescription({ "Store subsolution of non-linear solvers." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(settings, ECFMetaData()
			.setdescription({ "Store settings." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(FETI_data, ECFMetaData()
			.setdescription({ "Store FETI related data." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(catalyst, ECFMetaData()
			.setdescription({ "In-situ visualization by Catalyst." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(catalyst_sleep_time, ECFMetaData()
			.setdescription({ "The sleep time between each time steps when catalyst is used." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	collected = separate_bodies = separate_materials = false;
	REGISTER(collected, ECFMetaData()
			.setdescription({ "Collect the results to one file." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(separate_bodies, ECFMetaData()
			.setdescription({ "Separate bodies to regions." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(separate_materials, ECFMetaData()
			.setdescription({ "Separate materials to regions." })
			.setdatatype({ ECFDataType::BOOL }));

	domain_shrink_ratio = .95;
	cluster_shrink_ratio = .9;
	REGISTER(domain_shrink_ratio, ECFMetaData()
			.setdescription({ "Domain shrink ratio." })
			.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(cluster_shrink_ratio, ECFMetaData()
			.setdescription({ "Cluster shrink ratio." })
			.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(monitoring, ECFMetaData()
			.setdescription({ "Monitored region.", "STAT={AVG, MIN, MAX} PROPERTY"})
			.setdatatype({ ECFDataType::STRING, ECFDataType::STRING })
			.setpattern({ "REGION", "AVG TEMPERATURE" }));
}



