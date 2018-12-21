
#ifndef SRC_INPUT_WORKBENCH_COMMAND_H_
#define SRC_INPUT_WORKBENCH_COMMAND_H_

#include "parser.h"

namespace espreso {

struct Point;

struct NBlock: public WorkbenchParser {
	static size_t size;
	static const char* upper;
	static const char* lower;

	esint NUMFIELD, Solkey, NDMAX, NDSEL;

	esint lineSize, lineEndSize;
	esint indexSize, indexLength, valueSize, valueLength;

	NBlock();
	NBlock& parse(const char* begin);

	bool readData(std::vector<esint> &nIDs, std::vector<Point> &coordinates, double scaleFactor);

protected:
	bool index_x_y_z(std::vector<esint> &nIDs, std::vector<Point> &coordinates, double scaleFactor);
	bool index_solid_line_x_y_z(std::vector<esint> &nIDs, std::vector<Point> &coordinates, double scaleFactor);
};

}


#endif /* SRC_INPUT_WORKBENCH_COMMAND_H_ */
