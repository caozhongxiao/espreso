
#include "esinfo/systeminfo.h"
#include "basis/logging/logging.h"

#include <csignal>

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGTERM:
		ESINFO(espreso::ERROR) << "SIGTERM -termination request, sent to the program";
		break;
	case SIGSEGV:
		ESINFO(espreso::ERROR) << "SIGSEGV - invalid memory access (segmentation fault)";
		break;
	case SIGINT:
		ESINFO(espreso::ERROR) << "SIGINT - external interrupt, usually initiated by the user";
		break;
	case SIGILL:
		ESINFO(espreso::ERROR) << "SIGILL - invalid program image, such as invalid instruction";
		break;
	case SIGABRT:
		ESINFO(espreso::ERROR) << "SIGABRT - abnormal termination condition, as is e.g. initiated by std::abort()";
		break;
	case SIGFPE:
		ESINFO(espreso::ERROR) << "SIGFPE - erroneous arithmetic operation such as divide by zero";
		break;
	default:
		ESINFO(espreso::ERROR) << "ESPRESO trigger error " << signal << ".";
	}
}

namespace espreso {
namespace info {
namespace system {

void setSignals()
{
	switch (build) {
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

