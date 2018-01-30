
#include "../parser/esel.h"

#include "../../../basis/containers/tarray.h"
#include "../../../basis/utilities/parser.h"
#include "../../../basis/logging/logging.h"

using namespace espreso;

size_t ESel::size = 5;
const char* ESel::upper = "ESEL,";
const char* ESel::lower = "esel,";

ESel::ESel()
: type(Type::UNKNOWN), item(Item::UNKNOWN), comp(Comp::UNKNOWN),
  VMIN(0), VMAX(0), VINC(1),
  KABS(false)
{

}

ESel& ESel::parse(const char* begin)
{
	std::string commandLine = Parser::getLine(begin);
	std::vector<std::string> command = Parser::split(Parser::strip(commandLine), ",", false);

	switch (command.size()) {
	case 8:
		KABS = std::stoi(command[7]);
	case 7:
		if (command[6].size()) {
			VINC = std::stoi(command[6]);
		}
	case 6:
		if (command[5].size()) {
			size_t end;
			VMAX = std::stoi(command[5], &end);
			if (end != command[5].size()) {
				ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, VMAX='" << command[5] << "'";
			}
		}
	case 5:
		if (command[4].size()) {
			size_t end;
			VMIN = std::stoi(command[4], &end);
			if (end != command[4].size()) {
				ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, VMIN='" << command[4] << "'";
			}
		}
		if (command.size() < 6 || command[5].size() == 0) {
			VMAX = VMIN;
		}
	case 4:
		if (command[3].size()) {
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Comp='" << command[3] << "'";
		}
	case 3:
		if (StringCompare::caseInsensitiveEq("ELEM", command[2])) {
			item = Item::ELEM;
		}
		if (StringCompare::caseInsensitiveEq("ADJ", command[2])) {
			item = Item::ADJ;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("CENT", command[2])) {
			item = Item::CENT;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("TYPE", command[2])) {
			item = Item::TYPE;
		}
		if (StringCompare::caseInsensitiveEq("ENAME", command[2])) {
			item = Item::ENAME;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("MAT", command[2])) {
			item = Item::MAT;
		}
		if (StringCompare::caseInsensitiveEq("REAL", command[2])) {
			item = Item::REAL;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("ESYS", command[2])) {
			item = Item::ESYS;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("PART", command[2])) {
			item = Item::PART;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("LIVE", command[2])) {
			item = Item::LIVE;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("LAYER", command[2])) {
			item = Item::LAYER;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("SEC", command[2])) {
			item = Item::SEC;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("STRA", command[2])) {
			item = Item::STRA;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("SFE", command[2])) {
			item = Item::SFE;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("BFE", command[2])) {
			item = Item::BFE;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("PATH", command[2])) {
			item = Item::PATH;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
		if (StringCompare::caseInsensitiveEq("ETAB", command[2])) {
			item = Item::ETAB;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Item='" << command[2] << "'";
		}
	case 2:
		if (StringCompare::caseInsensitiveEq("S", command[1])) {
			type = Type::S;
		}
		if (StringCompare::caseInsensitiveEq("R", command[1])) {
			type = Type::R;
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Type='" << command[1] << "'";
		}
		if (StringCompare::caseInsensitiveEq("A", command[1])) {
			type = Type::A;
		}
		if (StringCompare::caseInsensitiveEq("U", command[1])) {
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Type='" << command[1] << "'";
			type = Type::U;
		}
		if (StringCompare::caseInsensitiveEq("ALL", command[1])) {
			type = Type::ALL;
		}
		if (StringCompare::caseInsensitiveEq("NONE", command[1])) {
			type = Type::NONE;
		}
		if (StringCompare::caseInsensitiveEq("INVE", command[1])) {
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Type='" << command[1] << "'";
			type = Type::INVE;
		}
		if (StringCompare::caseInsensitiveEq("STAT", command[1])) {
			ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel, Type='" << command[1] << "'";
			type = Type::STAT;
		}
		break;
	default:
		ESINFO(ERROR) << "ESPRESO Workbench parser error: unknown format of '" << commandLine << "'";
	}

	if (type == Type::UNKNOWN) {
		ESINFO(ERROR) << "ESPRESO Workbench parser error: unknown esel Type='" << command[1] << "'";
	}
	if (comp != Comp::UNKNOWN) {
		ESINFO(ERROR) << "ESPRESO Workbench parser error: implement esel Comp='" << command[3] << "'";
	}

	WorkbenchParser::fillIndices(begin, begin + commandLine.size());
	return *this;
}


