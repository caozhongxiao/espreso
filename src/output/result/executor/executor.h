
#ifndef SRC_OUTPUT_RESULT_EXECUTOR_EXECUTOR_H_
#define SRC_OUTPUT_RESULT_EXECUTOR_EXECUTOR_H_

#include "../resultstore.h"

#include <vector>

namespace espreso {

class ResultStoreExecutor: public ResultStoreBase {

public:
	virtual bool isCollected();
	virtual bool isSeparated();
	bool storeStep(const Step &step);

	virtual void addResultStore(ResultStoreBase *resultStore);
	virtual bool hasStore() { return _resultStore.size(); }

	ResultStoreExecutor(const Mesh &mesh, const OutputConfiguration &configuration): ResultStoreBase(mesh), _configuration(configuration) {}
	virtual ~ResultStoreExecutor();

protected:
	const OutputConfiguration &_configuration;

	std::vector<ResultStoreBase*> _resultStore;
};

}



#endif /* SRC_OUTPUT_RESULT_EXECUTOR_EXECUTOR_H_ */
