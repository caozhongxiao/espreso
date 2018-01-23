
#include "input.h"

#include "../../configuration.hpp"

using namespace espreso;

InputConfiguration::InputConfiguration()
{
	path = ".";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Path to an input example." })
			.setdatatype({ ECFDataType::STRING }));

	compress_numbers = false;
	REGISTER(compress_numbers, ECFMetaData()
			.setdescription({ "Make number compression (re-index nodes and element to make numbering continuous)." })
			.setdatatype({ ECFDataType::BOOL }));
}




