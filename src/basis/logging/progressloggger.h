
#ifndef SRC_BASIS_LOGGING_PROGRESSLOGGGER_H_
#define SRC_BASIS_LOGGING_PROGRESSLOGGGER_H_

#include "verbosity.h"
#include <cstdio>

namespace espreso {

class ProgressLogger: public Verbosity<ProgressLogger, 'v'> {

#define __ESLOG__COLOR__WHITE   "\x1b[0m"
#define __ESLOG__COLOR__RED     "\x1b[31m"
#define __ESLOG__COLOR__GREEN   "\x1b[32m"
#define __ESLOG__COLOR__YELLOW  "\x1b[33m"
#define __ESLOG__COLOR__BLUE    "\x1b[34m"
#define __ESLOG__COLOR__MAGENTA "\x1b[35m"
#define __ESLOG__COLOR__CYAN    "\x1b[36m"

public:
	void start(const char* region)
	{
		if (!rank) {
			printf("%*s%s", level, " ", region);
			fprintf(f, "%s", region);
		}
	}

	void checkpoint(const char* region)
	{
		if (!rank) {
			printf("%*s%s", level, " ", region);
			fprintf(f, "%s", region);
		}
	}

	void end(const char* region)
	{
		if (!rank) {
			printf("%*s%s", level, " ", region);
			fprintf(f, "%s", region);
		}
	}


	void param(const char* name, const int &value)
	{
		if (!rank) {
			printf(" [%s=%d]", name, value);
			fprintf(f, " [%s=%d]", name, value);
		}
	}

	void param(const char* name, const long &value)
	{
		if (!rank) {
			printf(" [%s=%ld]", name, value);
			fprintf(f, " [%s=%ld]", name, value);
		}
	}

	void param(const char* name, const long unsigned int &value)
	{
		if (!rank) {
			printf(" [%s=%ld]", name, value);
			fprintf(f, " [%s=%ld]", name, value);
		}
	}

	void param(const char* name, const double &value)
	{
		if (!rank) {
			printf(" [%s=%f]", name, value);
			fprintf(f, " [%s=%f]", name, value);
		}
	}

	void param(const char* name, const char* value)
	{
		if (!rank) {
			printf(" [%s=%s]", name, value);
			fprintf(f, " [%s=%s]", name, value);
		}
	}

	void ln()
	{
		if (!rank) {
			printf("\n");
			fprintf(f, "\n");
		}
	}

	void nextLoadStep(int step)
	{
		// do nothing
	}

	template <typename... Args>
	void info(const char* format, Args... args)
	{
		if (!rank) {
			printf(format, args...);
			fprintf(f, format, args...);
		}
	}

	template <typename... Args>
	void warning(const char* format, Args... args)
	{
		if (!rank) {
			printf(__ESLOG__COLOR__YELLOW);
			printf(format, args...);
			printf(__ESLOG__COLOR__WHITE);
			fprintf(f, format, args...);
		}
	}

	template <typename... Args>
	void debug(const char* format, Args... args)
	{
		if (!rank) {
			printf(__ESLOG__COLOR__BLUE);
			printf(format, args...);
			printf(__ESLOG__COLOR__WHITE);
			fprintf(f, format, args...);
		}
	}

	template <typename... Args>
	void error(const char* format, Args... args)
	{
		fprintf(stderr, __ESLOG__COLOR__RED);
		fprintf(stderr, format, args...);
		fprintf(stderr, __ESLOG__COLOR__WHITE);
		fprintf(f, format, args...);
		terminate();
	}

	ProgressLogger();
	~ProgressLogger();
	void setLogFile(const char* file);
	void terminate();

	int rank, size;
	FILE *f;
};

}


#endif /* SRC_BASIS_LOGGING_PROGRESSLOGGGER_H_ */
