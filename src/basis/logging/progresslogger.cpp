
#include "basis/logging/progresslogger.h"
#include "basis/utilities/sysutils.h"

#include <cstdlib>

using namespace espreso;

ProgressLogger::ProgressLogger()
: rank(0), size(1), f(NULL)
{

}

ProgressLogger::~ProgressLogger()
{
	if (f) {
		fclose(f);
	}
}

void ProgressLogger::setLogFile(const char* file)
{
	f = fopen(file, "wa");
}

void ProgressLogger::terminate()
{
	utils::printStack();
	fflush(stderr);
	fflush(f);
	exit(EXIT_FAILURE);
}


