
#ifndef SRC_INPUT_ABAQUS_PARSER_NSET_H_
#define SRC_INPUT_ABAQUS_PARSER_NSET_H_

#include "parser.h"

namespace espreso {

struct PlainAbaqusData;
struct Point;

struct Nset: public AbaqusParser {
	static size_t size;
	static const char* upper;
	static const char* lower;
	static const char* sentence;
	char NAME[MAX_NAME_SIZE];

	eslocal NUMFIELD, Solkey, NDMAX, NDSEL;

	eslocal lineSize, lineEndSize;
	eslocal indexSize, indexLength, valueSize, valueLength;

	Nset();
	Nset& parse(const char* begin);

	bool readData(std::vector<eslocal> &indices);


};

}


#endif /* SRC_INPUT_ABAQUS_COMMAND_H_ */
