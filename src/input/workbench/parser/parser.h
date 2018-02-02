
#ifndef SRC_INPUT_WORKBENCH_PARSER_PARSER_H_
#define SRC_INPUT_WORKBENCH_PARSER_PARSER_H_

#include <cstddef>
#include <vector>
#include <string>

#define MAX_NAME_SIZE 50 // upper bound on name size
#define MAX_COMMAND_SIZE 150 // upper bound on command size
#define MAX_LINE_SIZE 500 // upper bound on line size
#define MAX_LINE_STEP 3   // sometimes we need to red more lines to get full information

namespace espreso {

struct BlockEnd;

struct WorkbenchParser {
	static eslocal offset;
	static const char* begin;
	static const char* end;

	eslocal header;
	eslocal first, last;
	int fRank, lRank;

	WorkbenchParser()
	: header(-1),
	  first(-1), last(-1), fRank(-1), lRank(-1) {}

	void fillIndices(const char* header, const char* data);
	void fillIndices(const char* header, const char* first, const char* last);

	void fillDistribution(std::vector<BlockEnd> &blocksEnds, std::vector<eslocal> &distribution);
	const char* getFirst() const;
	const char* getLast() const;

	std::string command() const;

	void print(const char* data);
};
}



#endif /* SRC_INPUT_WORKBENCH_PARSER_PARSER_H_ */
