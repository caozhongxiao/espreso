
#ifndef SRC_CONFIG_ECF_PYTHONTESTGENERATOR_H_
#define SRC_CONFIG_ECF_PYTHONTESTGENERATOR_H_

#include "../configuration.h"

namespace espreso {

struct PythonTestGeneratorTable: public ECFObject {

	std::string columns, rows, value;

	PythonTestGeneratorTable();
};

struct PythonTestGenerator: public ECFObject {

	std::string output, env, run, exe;
	size_t measure_repetition, gather_level;

	std::map<size_t, std::string> levels, args;
	std::map<std::string, std::string> variables;

	std::map<std::string, PythonTestGeneratorTable> tables;

	PythonTestGenerator();
};

}


#endif /* SRC_CONFIG_ECF_PYTHONTESTGENERATOR_H_ */
