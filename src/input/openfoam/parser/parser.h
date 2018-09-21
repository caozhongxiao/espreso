
#ifndef SRC_INPUT_OPENFOAM_PARSER_PARSER_H_
#define SRC_INPUT_OPENFOAM_PARSER_PARSER_H_

#include <cstdlib>
#include <string>

namespace espreso {

struct OpenFOAMParser {

	const char *begin, *end, *current;
	char *endptr;

	OpenFOAMParser(const char *begin, const char *end);

	double readDouble()
	{
		double value = strtod(current, &endptr);
		current = endptr;
		return value;
	}

	double readInteger()
	{
		eslocal value = strtol(current, &endptr, 10);
		current = endptr;
		return value;
	}

	bool isEmpty()
	{
		return
				*current == ' ' ||
				*current == '\n' ||
				*current == '\r' ||
				*current == '\t';
	}

	std::string readString()
	{
		while (isEmpty()) { ++current; }
		const char *s = current;
		while (!isEmpty()) { ++current; }
		return std::string(s, current);
	}
};
}




#endif /* SRC_INPUT_OPENFOAM_PARSER_PARSER_H_ */
