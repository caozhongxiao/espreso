
#ifndef SRC_INPUT_WORKBENCH_PARSER_BLOCKEND_H_
#define SRC_INPUT_WORKBENCH_PARSER_BLOCKEND_H_

#include "../parser/parser.h"

namespace espreso {

struct BlockEnd: public WorkbenchParser {
	static size_t nSize;
	static size_t unixSize;
	static size_t winSize;

	static const char* nUpper;
	static const char* nLower;
	static const char* unixEnd;
	static const char* winEnd;

	BlockEnd();
	BlockEnd& parse(const char* begin);
};

}



#endif /* SRC_INPUT_WORKBENCH_PARSER_BLOCKEND_H_ */
