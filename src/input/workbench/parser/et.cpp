
#include "et.h"

#include "../../../basis/utilities/parser.h"
#include "../../../basis/logging/logging.h"

using namespace espreso;

size_t ET::size = 3;
const char* ET::upper = "ET,";
const char* ET::lower = "et,";

ET::ET()
: id(-1), type(-1)
{

}

ET& ET::parse(const char* begin)
{
	std::string commandLine = Parser::getLine(begin);
	std::vector<std::string> command = Parser::split(Parser::strip(commandLine), ",", false);

	if (command.size() < 3) {
		ESINFO(ERROR) << "ESPRESO Workbench parser error: unknown et format: " << commandLine;
	}

	if (
			!StringCompare::caseInsensitiveEq(command[1], "tid") &&
			!StringCompare::caseInsensitiveEq(command[1], "_tid") &&
			!StringCompare::caseInsensitiveEq(command[1], "cid")) {
		id = std::stoi(command[1]) - 1;
	}
	type = std::stoi(command[2]);

	WorkbenchParser::fillIndices(begin, begin);

	if (id != -1 && etype() == ETYPE::UNKNOWN) {
		ESINFO(ERROR) << "ESPRESO Workbench parser error: unknown et type: " << command[2];
	}

	return *this;
}


