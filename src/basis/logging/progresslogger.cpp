
#include "progressloggger.h"

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


