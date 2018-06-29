
#ifndef SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_
#define SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_

#include "../parser/parser.h"

namespace espreso {

struct PlainMeshData;
struct ET;

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

	void fixOffsets(std::vector<size_t> &dataOffsets);
	bool readData(const std::vector<ET> &et, PlainMeshData &mesh);

protected:
	bool solid(const std::vector<ET> &et, PlainMeshData &mesh);
	bool boundary(const std::vector<ET> &et, PlainMeshData &mesh);
};

}


#endif /* SRC_INPUT_WORKBENCH_PARSER_EBLOCK_H_ */
