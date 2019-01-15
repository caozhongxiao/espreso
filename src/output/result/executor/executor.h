
#ifndef SRC_OUTPUT_RESULT_EXECUTOR_EXECUTOR_H_
#define SRC_OUTPUT_RESULT_EXECUTOR_EXECUTOR_H_

#include "output/result/resultstore.h"

#include <vector>

namespace espreso {

class ResultStoreExecutor: public ResultStoreBase {

public:
	virtual bool isCollected();
	virtual bool isSeparated();
	bool storeStep();

	virtual void addResultStore(ResultStoreBase *resultStore);
	virtual bool hasStore() { return _resultStore.size(); }

	ResultStoreExecutor(const Mesh &mesh): ResultStoreBase(mesh) {}
	virtual ~ResultStoreExecutor();

protected:
	std::vector<ResultStoreBase*> _resultStore;
};

}



#endif /* SRC_OUTPUT_RESULT_EXECUTOR_EXECUTOR_H_ */
