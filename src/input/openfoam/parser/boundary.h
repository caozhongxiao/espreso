
#ifndef SRC_INPUT_OPENFOAM_PARSER_BOUNDARY_H_
#define SRC_INPUT_OPENFOAM_PARSER_BOUNDARY_H_

#include <string>
#include <vector>
#include <map>

#include "parser.h"

namespace espreso {

struct PlainOpenFOAMData;

struct OpenFOAMBoundary: public OpenFOAMParser {

	OpenFOAMBoundary(const char *begin, const char *end): OpenFOAMParser(begin, end) {}

	bool readData(std::map<std::string, std::vector<eslocal> > &eregions);
};

}



#endif /* SRC_INPUT_OPENFOAM_PARSER_BOUNDARY_H_ */
