
#include "parser.h"

using namespace espreso;

OpenFOAMParser::OpenFOAMParser(const char* begin, const char* end)
: begin(begin), end(end), current(begin), endptr(0)
{
	// TODO: parallel version

	while (*this->begin++ != '(');// { ++this->begin; }
//	if (*(this->begin + 1) == '\n' || *(this->begin + 1) == '\r') {
//		++this->begin;
//	}

	while (*this->end != ')') { --this->end; }
//	if (*(this->end + 1) == '\n' || *(this->end + 1) == '\r') {
//		--this->end;
//	}
}





