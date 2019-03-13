
#ifndef SRC_ESINFO_ESLOG_HPP_
#define SRC_ESINFO_ESLOG_HPP_

#include "eslog.h"
#include "basis/logging/progresslogger.h"

namespace espreso {
namespace eslog {

ProgressLogger& progress();

template <typename... Args>
void info(const char* format, Args... args)
{
	progress().info(format, args...);
}

template <typename... Args>
void solver(const char* format, Args... args)
{
	progress().info(format, args...);
}

template <typename... Args>
void linearsolver(const char* format, Args... args)
{
	if (progress().verbosity > 2) {
		progress().info(format, args...);
	}
}

template <typename... Args>
void duration(const char* format, Args... args)
{
	if (printtime()) {
		progress().info(format, args...);
	}
}

template <typename... Args>
void warning(const char* format, Args... args)
{
	progress().warning(format, args...);
}

template <typename... Args>
void debug(const char* format, Args... args)
{
	progress().debug(format, args...);
}

template <typename... Args>
void error(const char* format, Args... args)
{
	progress().error(format, args...);
}

template <typename... Args>
void globalerror(const char* format, Args... args)
{
	if (progress().rank == 0) {
		error(format, args...);
	}
}

}
}



#endif /* SRC_ESINFO_ESLOG_HPP_ */
