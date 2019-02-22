
#ifndef SRC_BASIS_LOGGING_OLDTIMELOGGER_H_
#define SRC_BASIS_LOGGING_OLDTIMELOGGER_H_

// TODO: dummy only for assign verebosity and m parameter

#include "verbosity.h"

namespace espreso {

class OldTimeLogger: public Verbosity<OldTimeLogger, 'm'> {

public:
	void start(const char* region)
	{

	}

	void checkpoint(const char* region)
	{

	}

	void end(const char* region)
	{

	}

	void param(const char* name, const int &value)
	{

	}

	void param(const char* name, const long &value)
	{

	}

	void param(const char* name, const long unsigned int &value)
	{

	}

	void param(const char* name, const double &value)
	{

	}

	void param(const char* name, const char* value)
	{
		// do nothing
	}

	void ln()
	{
		// do nothing
	}
};

}



#endif /* SRC_BASIS_LOGGING_OLDTIMELOGGER_H_ */
