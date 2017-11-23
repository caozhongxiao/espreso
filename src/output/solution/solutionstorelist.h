
#ifndef SRC_OUTPUT_SOLUTION_SOLUTIONSTORELIST_H_
#define SRC_OUTPUT_SOLUTION_SOLUTIONSTORELIST_H_

#include <vector>

#include "solutionstore.h"

namespace async { class Dispatcher; }

namespace espreso {

class OutputConfiguration;

class SolutionStoreList: public Store {

public:
	static SolutionStoreList* createAsynchronizedStore(const Mesh &mesh, const OutputConfiguration &configuration);
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

	~SolutionStoreList()
	{
		for (size_t i = 0; i < _results.size(); ++i) {
			delete _results[i];
		}
	}

protected:
	void add(SolutionStore* store) { _results.push_back(store); }

	std::vector<SolutionStore*> _results;

	static SolutionStoreList *_resultStoreList;
	static async::Dispatcher *_dispatcher;
};

}



#endif /* SRC_OUTPUT_SOLUTION_SOLUTIONSTORELIST_H_ */
