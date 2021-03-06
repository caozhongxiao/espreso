
#ifndef SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_
#define SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_

#include "parser.h"

namespace espreso {

struct PlainWorkbenchData;
struct ET;

struct EBlock: public WorkbenchParser {
	static size_t size;
	static const char* upper;
	static const char* lower;

	esint NUM_NODES, Solkey, NDMAX, NDSEL;

	esint lineSize, elementSize, lineEndSize;
	esint valueSize, valueLength;

	EBlock();
	EBlock& parse(const char* begin);

	bool solidBlock();

	void fixOffsets(std::vector<size_t> &dataOffsets);
	bool readData(const std::vector<ET> &et, PlainWorkbenchData &mesh);

protected:
	bool solid(const std::vector<ET> &et, PlainWorkbenchData &mesh);
	bool boundary(const std::vector<ET> &et, PlainWorkbenchData &mesh);
};

}


#endif /* SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_ */
