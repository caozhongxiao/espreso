
#ifndef SRC_INPUT_OPENFOAM_PARSER_BOUNDARY_H_
#define SRC_INPUT_OPENFOAM_PARSER_BOUNDARY_H_

#include "parser.h"

#include <vector>

namespace espreso {

#define MAX_NAME_SIZE 80

struct OpenFOAMBoundaryData {
	char name[MAX_NAME_SIZE];
	size_t startFace, nFaces;

	OpenFOAMBoundaryData();
	OpenFOAMBoundaryData(const std::string &name, size_t startFace, size_t nFaces);
};

struct OpenFOAMBoundary: public OpenFOAMSeparateParser {

	OpenFOAMBoundary(const char *begin, const char *end): OpenFOAMSeparateParser(begin, end) {}

	bool readData(std::vector<OpenFOAMBoundaryData> &boundaries);
};

}



#endif /* SRC_INPUT_OPENFOAM_PARSER_BOUNDARY_H_ */
