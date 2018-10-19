
#ifndef SRC_INPUT_OPENFOAM_PARSER_FACES_H_
#define SRC_INPUT_OPENFOAM_PARSER_FACES_H_

#include <vector>

#include "parser.h"

namespace espreso {

struct PlainOpenFOAMData;

struct OpenFOAMFaces: public OpenFOAMCollectiveParser {

	OpenFOAMFaces(const char *begin, const char *end): OpenFOAMCollectiveParser(begin, end) {}

	bool readFaces(PlainOpenFOAMData &data);
	bool readParents(std::vector<eslocal> &data);
};

}



#endif /* SRC_INPUT_OPENFOAM_PARSER_FACES_H_ */
