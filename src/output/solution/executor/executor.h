
#ifndef SRC_OUTPUT_SOLUTION_EXECUTOR_EXECUTOR_H_
#define SRC_OUTPUT_SOLUTION_EXECUTOR_EXECUTOR_H_

#include "../solutionstore.h"

namespace espreso {

class SolutionStoreExecutor: public SolutionStore {

public:
	void storePreprocessedData();
	bool storeSolution(const Step &step);

	SolutionStoreExecutor(const Mesh &mesh, const OutputConfiguration &configuration): SolutionStore(mesh, configuration) {}
};

}



#endif /* SRC_OUTPUT_SOLUTION_EXECUTOR_EXECUTOR_H_ */
