
#ifndef SRC_OUTPUT_RESULTS_EXECUTOR_DIRECTEXECUTOR_H_
#define SRC_OUTPUT_RESULTS_EXECUTOR_DIRECTEXECUTOR_H_

#include "executor.h"

namespace espreso {

class DirectExecutor: public ResultStoreExecutor {

public:
	void updateMesh();
	void updateSolution();

	DirectExecutor(const OutputConfiguration &configuration): ResultStoreExecutor(configuration) {}
};

}



#endif /* SRC_OUTPUT_RESULTS_EXECUTOR_DIRECTEXECUTOR_H_ */
