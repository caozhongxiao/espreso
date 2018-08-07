
#include "../parser/blockend.h"

using namespace espreso;

size_t BlockEnd::nSize = 2;
const char* BlockEnd::nUpper = "N,";
const char* BlockEnd::nLower = "n,";

size_t BlockEnd::unixSize = 3;
const char* BlockEnd::unixEnd = "-1\n";

size_t BlockEnd::winSize = 4;
const char* BlockEnd::winEnd = "-1\r\n";

BlockEnd::BlockEnd()
{

}

BlockEnd& BlockEnd::parse(const char* begin)
{
	while (*(--begin) != '\n');
	WorkbenchParser::fillIndices(begin + 1, begin + 1);
	return *this;
}

