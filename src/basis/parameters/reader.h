
#ifndef SRC_BASIS_PARAMETERS_READER_H_
#define SRC_BASIS_PARAMETERS_READER_H_

#include "mpi.h"

#include <getopt.h>

#include "esconfig.h"
#include "parameter.h"
#include "../logging/logging.h"


namespace espreso {

class ParametersReader {

public:
	static Configuration arguments(int *argc, char*** argv, const std::vector<Parameter> &params = config::parameters);
	static Configuration configuration(const Configuration &conf, const std::vector<Parameter> &params = config::parameters);
	static Configuration pickConfiguration(const Configuration &conf, const std::vector<Parameter> &params = config::parameters);

	static void printParameters(const std::vector<Parameter> &params, size_t verboseLevel);
	static void printParametersHelp(const std::vector<Parameter> &params, size_t verboseLevel);

protected:
	ParametersReader(const std::vector<Parameter> &parameters);
	Configuration read(const Configuration &configuration, size_t verboseLevel);

	std::vector<Parameter> _parameters;

private:
	static void printHelp(size_t verboseLevel);
	bool setParameter(const std::string &parameter, const std::string &value);
};

}



#endif /* SRC_BASIS_PARAMETERS_READER_H_ */