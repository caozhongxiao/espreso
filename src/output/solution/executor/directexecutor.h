
#ifndef SRC_OUTPUT_SOLUTION_EXECUTOR_DIRECTEXECUTOR_H_
#define SRC_OUTPUT_SOLUTION_EXECUTOR_DIRECTEXECUTOR_H_

#include "executor.h"

namespace espreso {

class DirectExecutor: public SolutionStoreExecutor {

public:
	void updateMesh();
	void updateSolution();

	DirectExecutor(const Mesh &mesh, const OutputConfiguration &configuration): SolutionStoreExecutor(mesh, configuration) {}
};

}



#endif /* SRC_OUTPUT_SOLUTION_EXECUTOR_DIRECTEXECUTOR_H_ */
