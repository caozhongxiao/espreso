
#include "systeminfo.h"
#include "eslog.h"

#include <csignal>

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGTERM:
		espreso::eslog::error("SIGTERM -termination request, sent to the program.\n");
		break;
	case SIGSEGV:
		espreso::eslog::error("SIGSEGV - invalid memory access (segmentation fault).\n");
		break;
	case SIGINT:
		espreso::eslog::error("SIGINT - external interrupt, usually initiated by the user.\n");
		break;
	case SIGILL:
		espreso::eslog::error("SIGILL - invalid program image, such as invalid instruction.\n");
		break;
	case SIGABRT:
		espreso::eslog::error("SIGABRT - abnormal termination condition, as is e.g. initiated by std::abort().\n");
		break;
	case SIGFPE:
		espreso::eslog::error("SIGFPE - erroneous arithmetic operation such as divide by zero.\n");
		break;
	default:
		espreso::eslog::error("ESPRESO trigger an error.\n");
	}
}

namespace espreso {
namespace info {
namespace system {

OPERATIONSYSTEM os()
{
	return OPERATIONSYSTEM::UNIX;
}

INSTRUCTIONSET instructionSet()
{
	return INSTRUCTIONSET::SSE;
}

BUILD build()
{
#ifdef __ESMODE__
	return BUILD::__ESMODE__;
#else
	return BUILD::RELEASE;
#endif
}

const char* commit()
{
#ifdef __ESCOMMIT__
	return __ESCOMMIT__;
#else
	return "undefined";
#endif
}

void setSignals()
{
	switch (build()) {
	case BUILD::RELEASE:
	case BUILD::MEASUREMENT:
	case BUILD::DEVEL:
	case BUILD::DEBUG:
		std::signal(SIGTERM, signalHandler);
		std::signal(SIGSEGV, signalHandler);
		std::signal(SIGINT, signalHandler);
		std::signal(SIGILL, signalHandler);
		std::signal(SIGABRT, signalHandler);
		std::signal(SIGFPE, signalHandler);
		break;
	}
}

}
}
}

