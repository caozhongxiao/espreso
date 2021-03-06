
#ifndef SRC_BASIS_LOGGING_PROGRESSLOGGER_H_
#define SRC_BASIS_LOGGING_PROGRESSLOGGER_H_

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
	void start(const char* region, const char* section)
	{
		if (!grank) {
			printf("%*s%s", level, " ", region);
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, "%*s%s", level, " ", region);
			fflush(f);
		}
	}

	void checkpoint(const char* region)
	{
		if (!grank) {
			printf("%*s%s", level, " ", region);
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, "%*s%s", level, " ", region);
			fflush(f);
		}
	}

	void end(const char* region)
	{
		if (!grank) {
			printf("%*s%s", level, " ", region);
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, "%*s%s", level, " ", region);
			fflush(f);
		}
	}


	void param(const char* name, const int &value)
	{
		if (!grank) {
			printf(" [%s=%d]", name, value);
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, " [%s=%d]", name, value);
			fflush(f);
		}
	}

	void param(const char* name, const long &value)
	{
		if (!grank) {
			printf(" [%s=%ld]", name, value);
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, " [%s=%ld]", name, value);
			fflush(f);
		}
	}

	void param(const char* name, const long unsigned int &value)
	{
		if (!grank) {
			printf(" [%s=%lu]", name, value);
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, " [%s=%lu]", name, value);
			fflush(f);
		}
	}

	void param(const char* name, const double &value)
	{
		if (!grank) {
			printf(" [%s=%f]", name, value);
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, " [%s=%f]", name, value);
			fflush(f);
		}
	}

	void param(const char* name, const char* value)
	{
		if (!grank) {
			printf(" [%s=%s]", name, value);
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, " [%s=%s]", name, value);
			fflush(f);
		}
	}

	void ln()
	{
		if (!grank) {
			printf("\n");
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, "\n");
			fflush(f);
		}
	}

	void nextLoadStep(int step)
	{
		// do nothing
	}

	template <typename... Args>
	void info(const char* format, Args... args)
	{
		if (!grank) {
			printf(format, args...);
			fflush(stdout);
		}
		if (!rank) {
			fprintf(f, format, args...);
			fflush(f);
		}
	}

	template <typename... Args>
	void warning(const char* format, Args... args)
	{
		if (!grank) {
			char* colored = getColored(format, __ESLOG__COLOR__YELLOW);
			printf(colored, args...);
			fflush(stdout);
			delete[] colored;
		}
		if (!rank) {
			fprintf(f, format, args...);
			fflush(f);
		}
	}

	template <typename... Args>
	void debug(const char* format, Args... args)
	{
		if (!grank) {
			char* colored = getColored(format, __ESLOG__COLOR__BLUE);
			printf(colored, args...);
			fflush(stdout);
			delete[] colored;
		}
		if (!rank) {
			fprintf(f, format, args...);
			fflush(f);
		}
	}

	template <typename... Args>
	void error(const char* format, Args... args)
	{
		char* colored = getColored(format, __ESLOG__COLOR__RED);
		fprintf(stderr, colored, args...);
		fflush(stderr);
		delete[] colored;
		fprintf(f, format, args...);
		fflush(f);
		terminate();
	}

	char* getColored(const char* str, const char *color);

	ProgressLogger();
	~ProgressLogger();
	void setLogFile(const char* file);
	void terminate();

	int rank, size;
	int grank, gsize;
	FILE *f;
};

}


#endif /* SRC_BASIS_LOGGING_PROGRESSLOGGER_H_ */
