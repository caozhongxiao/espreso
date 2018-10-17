
#include "boundary.h"

#include "../../../basis/utilities/parser.h"

#include <cstring>
#include <numeric>

using namespace espreso;

OpenFOAMBoundaryData::OpenFOAMBoundaryData()
: startFace(0), nFaces(0)
{
	memset(name, '\0', MAX_NAME_SIZE);
}

OpenFOAMBoundaryData::OpenFOAMBoundaryData(const std::string &name, size_t startFace, size_t nFaces)
: startFace(startFace), nFaces(nFaces)
{
	memset(this->name, '\0', MAX_NAME_SIZE);
	memcpy(this->name, name.data(), name.size() < MAX_NAME_SIZE ? name.size() : MAX_NAME_SIZE);
}

bool OpenFOAMBoundary::readData(std::vector<OpenFOAMBoundaryData> &boundaries)
{
	const char *c = begin;

	size_t nFaces, startFace;
	std::string name, parameter;

	while (*c != ')') {
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

		boundaries.push_back({name, startFace, nFaces});
	}
	return true;
}


