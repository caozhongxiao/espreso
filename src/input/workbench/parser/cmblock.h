
#ifndef SRC_INPUT_WORKBENCH_PARSER_CMBLOCK_H_
#define SRC_INPUT_WORKBENCH_PARSER_CMBLOCK_H_

#include "../parser/parser.h"

namespace espreso {

struct CMBlock: public WorkbenchParser {
	enum class Entity: int {
		NODE,
		ELEMENT
	};

	static size_t size;
	static const char* upper;
	static const char* lower;

	char name[21];
	Entity entity;
	eslocal NUMITEMS;

	eslocal lineSize, lineEndSize;
	eslocal valueSize, valueLength;

	CMBlock();
	CMBlock& parse(const char* begin);

	bool readData(std::vector<eslocal> &indices);
};

}


#endif /* SRC_INPUT_WORKBENCH_PARSER_CMBLOCK_H_ */
