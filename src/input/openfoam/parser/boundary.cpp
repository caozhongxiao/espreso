
#include "boundary.h"
#include "../openfoam.h"

#include "../../../basis/utilities/parser.h"

#include <numeric>

using namespace espreso;

bool OpenFOAMBoundary::readData(std::map<std::string, std::vector<eslocal> > &eregions)
{
	current = begin;

	size_t nFaces, startFace;
	std::string name, parameter;

	while (*current != ')') {
		name = readString();

		while (*current++ != '{');
		while (*current != '}') {
			parameter = readString();
			if (StringCompare::caseSensitiveEq(parameter, "nFaces")) {
				nFaces = readInteger();
			}
			if (StringCompare::caseSensitiveEq(parameter, "startFace")) {
				startFace = readInteger();
			}
			while (*current == ';' || isEmpty()) { ++current; }
		}
		++current;
		while (isEmpty()) { ++current; }

//		std::vector<eslocal> &elements = eregions[name];
//		elements.resize(nFaces);
//		std::iota(elements.begin(), elements.end(), startFace);
	}
	return true;
}


