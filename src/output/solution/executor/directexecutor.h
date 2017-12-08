
#ifndef SRC_OUTPUT_SOLUTION_EXECUTOR_DIRECTEXECUTOR_H_
#define SRC_OUTPUT_SOLUTION_EXECUTOR_DIRECTEXECUTOR_H_

#include "executor.h"

namespace espreso {

class EnSight;

class DirectExecutor: public SolutionStoreExecutor {

public:
	void updateMesh();
	void updateSolution(const Step &step);

	DirectExecutor(const Mesh &mesh, const OutputConfiguration &configuration);
	~DirectExecutor();

protected:
	EnSight *_ensight;
};

}



#endif /* SRC_OUTPUT_SOLUTION_EXECUTOR_DIRECTEXECUTOR_H_ */
