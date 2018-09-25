
#ifndef SRC_CONFIG_ECF_PROCESSESREDUCTION_H_
#define SRC_CONFIG_ECF_PROCESSESREDUCTION_H_

#include <cstddef>

namespace espreso {

struct ProcessesReduction {
	enum class Granularity {
		NODES,
		PROCESSES
	};

	enum class Pattern {
		PREFIX,
		SUBSET
	};

	Granularity granularity= Granularity::PROCESSES;
	Pattern pattern = Pattern::SUBSET;
	size_t reduction_ratio = 1;
};

}



#endif /* SRC_CONFIG_ECF_PROCESSESREDUCTION_H_ */
