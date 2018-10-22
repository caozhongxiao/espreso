
#ifndef SRC_INPUT_OPENFOAM_PARSER_ZONES_H_
#define SRC_INPUT_OPENFOAM_PARSER_ZONES_H_

#include <vector>

#include "parser.h"

namespace espreso {

struct PlainOpenFOAMData;
struct ParallelFile;

struct OpenFOAMZones: public OpenFOAMCollectiveParser {

	OpenFOAMZones(ParallelFile &pfile);

	int getZones();
	void synchronize(int zones, std::vector<char> &names, std::vector<size_t> &offsets);
	void readData(std::vector<eslocal> &indices, size_t begin, size_t end);

	bool readPoints(PlainOpenFOAMData &data);
	bool readFaces(PlainOpenFOAMData &data);
	bool readCells(PlainOpenFOAMData &data);

	ParallelFile &_pfile;
};

}


#endif /* SRC_INPUT_OPENFOAM_PARSER_ZONES_H_ */
