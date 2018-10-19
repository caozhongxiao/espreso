
#include "parser.h"

#include "mpi.h"
#include "../../../config/ecf/environment.h"
#include "../../../basis/utilities/communication.h"

using namespace espreso;

OpenFOAMParser::OpenFOAMParser(const char* begin, const char* end)
: begin(begin), end(end)
{

}

OpenFOAMSeparateParser::OpenFOAMSeparateParser(const char *begin, const char *end)
: OpenFOAMParser(begin, end)
{
	while (*this->begin++ != '(');
	while (*this->end-- != ')');
}

OpenFOAMCollectiveParser::OpenFOAMCollectiveParser(const char *begin, const char *end)
: OpenFOAMParser(begin, end)
{
	int found, min, max;
	const char *p;

	p = begin;
	while (*p != '(' && p < end) { ++p; };
	if (p == end) {
		found = environment->MPIsize;
	} else {
		found = environment->MPIrank;
	}

	MPI_Allreduce(&found, &min, 1, MPI_INT, MPI_MIN, environment->MPICommunicator);
	if (environment->MPIrank == min) {
		this->begin = p + 1;
	}
	if (environment->MPIrank < min) {
		this->begin = p;
	}

	p = end;
	while (*p != ')' && begin < p) { --p; }
	if (p == begin) {
		found = 0;
	} else {
		found = environment->MPIrank;
	}

	MPI_Allreduce(&found, &max, 1, MPI_INT, MPI_MAX, environment->MPICommunicator);
	if (max <= environment->MPIrank) {
		this->end = p;
	}
}





