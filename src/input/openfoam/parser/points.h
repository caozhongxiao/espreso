
#ifndef SRC_INPUT_OPENFOAM_PARSER_POINTS_H_
#define SRC_INPUT_OPENFOAM_PARSER_POINTS_H_

#include <vector>

#include "parser.h"

namespace espreso {

struct Point;

struct OpenFOAMPoints: public OpenFOAMParser {

	OpenFOAMPoints(const char *begin, const char *end): OpenFOAMParser(begin, end) {}

	bool readData(std::vector<eslocal> &nIDs, std::vector<Point> &coordinates, double scaleFactor);
};

}


#endif /* SRC_INPUT_OPENFOAM_PARSER_POINTS_H_ */
