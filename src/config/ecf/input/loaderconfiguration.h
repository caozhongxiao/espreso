
#ifndef SRC_CONFIG_ECF_INPUT_LOADERCONFIGURATION_H_
#define SRC_CONFIG_ECF_INPUT_LOADERCONFIGURATION_H_

#include <cstddef>

namespace espreso {

struct LoaderConfiguration {
	enum class LoadBy {
		NODES,
		PROCESSES
	};

	enum class LoadPattern {
		PREFIX,
		SUBSET
	};


	LoadBy load_by = LoadBy::PROCESSES;
	LoadPattern load_pattern = LoadPattern::SUBSET;
	size_t pattern_value = 1;
};

}



#endif /* SRC_CONFIG_ECF_INPUT_LOADERCONFIGURATION_H_ */
