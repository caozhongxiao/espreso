
#include "basis/logging/progresslogger.h"
#include "basis/utilities/sysutils.h"

#include <cstdlib>
#include <cstring>

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

char* ProgressLogger::getColored(const char* str, const char *color)
{
	char *colored = new char[strlen(str) + 15];
	strcpy(colored, color);
	strcat(colored, str);
	strcat(colored, __ESLOG__COLOR__WHITE);
	return colored;
}

void ProgressLogger::terminate()
{
	utils::printStack();
	fflush(stderr);
	fflush(f);
	exit(EXIT_FAILURE);
}


