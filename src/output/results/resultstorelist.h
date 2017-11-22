
#ifndef SRC_OUTPUT_RESULTS_RESULTSTORELIST_H_
#define SRC_OUTPUT_RESULTS_RESULTSTORELIST_H_

#include <vector>

#include "resultstore.h"

namespace async { class Dispatcher; }

namespace espreso {

class OutputConfiguration;

class ResultStoreList: public Store {

public:
	static ResultStoreList* createAsynchronizedStore(const OutputConfiguration &configuration);
	static void destroyAsynchronizedStore();
	static bool isStoreNode();
	static bool isComputeNode();

	void storePreprocessedData()
	{
		for (size_t i = 0; i < _results.size(); ++i) {
			_results[i]->storePreprocessedData();
		}
	}

	void updateMesh()
	{
		for (size_t i = 0; i < _results.size(); ++i) {
			_results[i]->updateMesh();
		}
	}

	void updateSolution()
	{
		for (size_t i = 0; i < _results.size(); ++i) {
			_results[i]->updateSolution();
		}
	}

	~ResultStoreList()
	{
		for (size_t i = 0; i < _results.size(); ++i) {
			delete _results[i];
		}
	}

protected:
	void add(ResultStore* store) { _results.push_back(store); }

	std::vector<ResultStore*> _results;

	static ResultStoreList *_resultStoreList;
	static async::Dispatcher *_dispatcher;
};

}



#endif /* SRC_OUTPUT_RESULTS_RESULTSTORELIST_H_ */
