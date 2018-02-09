
#ifndef SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_
#define SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_

#include "../parser/parser.h"

namespace espreso {

struct EData;

struct EBlock: public WorkbenchParser {
	static size_t size;
	static const char* upper;
	static const char* lower;

	eslocal NUM_NODES, Solkey, NDMAX, NDSEL;

	eslocal lineSize, elementSize, lineEndSize;
	eslocal valueSize, valueLength;

	EBlock();
	EBlock& parse(const char* begin);

	bool solidBlock();

	void fixOffsets(std::vector<eslocal> &dataOffsets);
	bool readSolid(std::vector<eslocal> &esize, std::vector<eslocal> &nodes, std::vector<EData> &data);
	bool readBoundary(std::vector<eslocal> &esize, std::vector<eslocal> &nodes, std::vector<EData> &edata);

protected:
	bool solid(std::vector<eslocal> &esize, std::vector<eslocal> &nodes, std::vector<EData> &data);
	bool boundary(std::vector<eslocal> &esize, std::vector<eslocal> &nodes, std::vector<EData> &edata);
};

}


#endif /* SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_ */
