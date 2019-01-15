
#include "decomposition.h"
#include "config/configuration.hpp"

espreso::METISConfiguration::METISConfiguration()
{
	objective_type = OBJECTIVE_TYPE::VOLUME;
	REGISTER(objective_type, ECFMetaData()
			.setdescription({ "Decomposition respect materials." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("VOLUME").setdescription("METIS tries to minimize communication volume. Recommended to examples with various element size."))
			.addoption(ECFOption().setname("EDGECUT").setdescription("METIS tries to minimize edgecut.")));


	refinement = false;
	REGISTER(refinement, ECFMetaData()
			.setdescription({ "Apply refinement." })
			.setdatatype({ ECFDataType::BOOL }));

	ecfdescription->addSeparator();

	granularity = ProcessesReduction::Granularity::PROCESSES;
	REGISTER(granularity, ECFMetaData()
			.setdescription({ "A decomposer granularity. " })
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

espreso::DecompositionConfiguration::DecompositionConfiguration()
{
	path = "EBF";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Path to store preprocessed mesh in ESPRESO binary format (*.ebf)." })
			.setdatatype({ ECFDataType::STRING }));

	ecfdescription->addSpace();

	mpi_procs = 0;
	domains = 0;
	REGISTER(mpi_procs, ECFMetaData()
			.setdescription({ "Number of MPI processes of *.ebf data." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
	REGISTER(domains, ECFMetaData()
			.setdescription({ "Number of domains for each cluster (Keep 0 for automatic decomposition)." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	ecfdescription->addSpace();

	balance_clusters = false;
	REGISTER(balance_clusters, ECFMetaData()
			.setdescription({ "Balance elements among MPI processes." })
			.setdatatype({ ECFDataType::BOOL }));

	ecfdescription->addSpace();

	separate_materials = false;
	separate_regions = false;
	separate_etypes = false;
	REGISTER(separate_materials, ECFMetaData()
			.setdescription({ "Decomposition respect materials." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(separate_regions, ECFMetaData()
			.setdescription({ "Decomposition respect regions." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(separate_etypes, ECFMetaData()
			.setdescription({ "Decomposition respect elements types." })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(metis_options, ECFMetaData()
			.setdescription({ "ParMETIS/METIS options." }));
}




