
#include "progressloggger.h"
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
	std::vector<std::string> stack;
	utils::currentStack(stack);
	for (size_t i = 0; i < stack.size(); i++) {
		info("%s", stack[i].c_str());
	}
	fflush(stderr);
	fflush(f);
	exit(EXIT_FAILURE);
}


