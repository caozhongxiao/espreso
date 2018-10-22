
#include "boundary.h"

#include "../openfoam.h"

#include "../../../basis/utilities/parser.h"

using namespace espreso;

bool OpenFOAMBoundary::readData(PlainOpenFOAMData &data)
{
	size_t nFaces, startFace;
	std::string name, parameter;

	const char *c = begin - 3;
	while (*c != '\n') { c--; } // go before number of boundaries

	int n = readInteger(c);
	while (*c++ != '('); // skip '('

	for (int i = 0; i < n; i++) {
		name = readString(c);

		while (*c++ != '{');
		while (*c != '}') {
			parameter = readString(c);
			if (StringCompare::caseSensitiveEq(parameter, "nFaces")) {
				nFaces = readInteger(c);
			}
			if (StringCompare::caseSensitiveEq(parameter, "startFace")) {
				startFace = readInteger(c);
			}
			while (*c == ';' || isEmpty(c)) { ++c; }
		}
		++c;
		while (isEmpty(c)) { ++c; }

		auto &indices = data.eregions[data.boundaryprefix + name];
		for (size_t f = 0; f < data.fIDs.size(); f++) {
			if (startFace <= data.fIDs[f] && data.fIDs[f] < startFace + nFaces) {
				indices.push_back(data.fIDs[f]);
			}
		}
	}

	return true;
}


