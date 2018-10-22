
#ifndef SRC_INPUT_OPENFOAM_PARSER_BOUNDARY_H_
#define SRC_INPUT_OPENFOAM_PARSER_BOUNDARY_H_

#include "parser.h"

#include <vector>

namespace espreso {

struct PlainOpenFOAMData;

struct OpenFOAMBoundary: public OpenFOAMSeparateParser {

	OpenFOAMBoundary(const char *begin, const char *end): OpenFOAMSeparateParser(begin, end) {}

	bool readData(PlainOpenFOAMData &data);
};

}



#endif /* SRC_INPUT_OPENFOAM_PARSER_BOUNDARY_H_ */
