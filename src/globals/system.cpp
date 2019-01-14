
#include "system.h"

#include "basis/logging/logging.h"

#include <csignal>

using namespace espreso;

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGTERM:
		ESINFO(ERROR) << "SIGTERM -termination request, sent to the program";
		break;
	case SIGSEGV:
		ESINFO(ERROR) << "SIGSEGV - invalid memory access (segmentation fault)";
		break;
	case SIGINT:
		ESINFO(ERROR) << "SIGINT - external interrupt, usually initiated by the user";
		break;
	case SIGILL:
		ESINFO(ERROR) << "SIGILL - invalid program image, such as invalid instruction";
		break;
	case SIGABRT:
		ESINFO(ERROR) << "SIGABRT - abnormal termination condition, as is e.g. initiated by std::abort()";
		break;
	case SIGFPE:
		ESINFO(ERROR) << "SIGFPE - erroneous arithmetic operation such as divide by zero";
		break;
	default:
		ESINFO(ERROR) << "ESPRESO trigger error " << signal << ".";
	}
}

void system::setSignals()
{
	switch (build()) {
	case BUILD::DEBUG:
	case BUILD::RELEASE:
		std::signal(SIGTERM, signalHandler);
		std::signal(SIGSEGV, signalHandler);
		std::signal(SIGINT, signalHandler);
		std::signal(SIGILL, signalHandler);
		std::signal(SIGABRT, signalHandler);
		std::signal(SIGFPE, signalHandler);
		break;
	}
}
