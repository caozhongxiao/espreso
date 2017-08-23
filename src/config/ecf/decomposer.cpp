
#include "decomposer.h"
#include "../configuration.hpp"

espreso::DecomposerConfiguration::DecomposerConfiguration()
{
	parts = "4";
	prefix = "ESDATA";

	REGISTER(parts, ECFMetaData()
			.setdescription({ "Each MPI process will be decomposed into the specified number of parts (e.q. 1 2 4)." })
			.setdatatype({ ECFDataType::STRING }));
	REGISTER(prefix, ECFMetaData()
			.setdescription({ "Decomposition will be saved into PREFIX{PARTS} directories." })
			.setdatatype({ ECFDataType::STRING }));
}



