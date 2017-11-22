
#ifndef SRC_OUTPUT_RESULTS_EXECUTOR_EXECUTOR_H_
#define SRC_OUTPUT_RESULTS_EXECUTOR_EXECUTOR_H_

#include "../resultstore.h"

namespace espreso {

class ResultStoreExecutor: public ResultStore {

public:
	void storePreprocessedData();

	ResultStoreExecutor(const OutputConfiguration &configuration): ResultStore(configuration) {}
};

}



#endif /* SRC_OUTPUT_RESULTS_EXECUTOR_EXECUTOR_H_ */
