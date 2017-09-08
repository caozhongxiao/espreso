
#include "input.h"

#include "../../configuration.hpp"

using namespace espreso;

InputConfiguration::InputConfiguration()
{
	path = ".";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Path to an input example." })
			.setdatatype({ ECFDataType::STRING }));

	domains = 4;
	REGISTER(domains, ECFMetaData()
			.setdescription({ "Number of domains for each cluster (MPI process)." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
}




