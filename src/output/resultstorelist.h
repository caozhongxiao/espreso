
#ifndef SRC_OUTPUT_RESULTSTORELIST_H_
#define SRC_OUTPUT_RESULTSTORELIST_H_

//namespace async { class Dispatcher; }

namespace espreso {

class OutputConfiguration;

class ResultStoreList {

public:
	static ResultStoreList* createAsynchronizedStore(const OutputConfiguration &configuration) { return NULL; }
	static void destroyAsynchronizedStore() {}
	static bool isStoreNode() { return false; }
	static bool isComputeNode() { return !isStoreNode(); }

	ResultStoreList() { }

	~ResultStoreList()
	{
//		for (size_t i = 0; i < _results.size(); ++i) {
//			delete _results[i];
//		}
	}

	void updateMesh()
	{
//		for (size_t i = 0; i < _results.size(); ++i) {
//			_results[i]->updateMesh();
//		}
	}

	void updateSolution()
	{
//		for (size_t i = 0; i < _results.size(); ++i) {
//			_results[i]->updateSolution();
//		}
	}

protected:
//	std::vector<ResultStore*> _results;
//
//	static ResultStoreList *_resultStoreList;
//	static async::Dispatcher *_dispatcher;
};

}



#endif /* SRC_OUTPUT_RESULTSTORELIST_H_ */
