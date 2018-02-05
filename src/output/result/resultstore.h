
#ifndef SRC_OUTPUT_RESULT_RESULTSTORE_H_
#define SRC_OUTPUT_RESULT_RESULTSTORE_H_

#include <string>

namespace async { class Dispatcher; }

namespace espreso {

struct Step;
class Mesh;
class OutputConfiguration;
class ResultStoreExecutor;

class ResultStoreBase {

public:
	virtual bool isCollected() =0;
	virtual bool isDistributed() =0;

	virtual void updateMesh() =0;
	virtual void updateSolution(const Step &step) =0;

	virtual const Mesh& mesh() const { return _mesh; }

	virtual ~ResultStoreBase() {};

protected:
	ResultStoreBase(const Mesh &mesh);

	const Mesh &_mesh;
	std::string _directory;
};

class ResultStore {

public:
	static ResultStore* createAsynchronizedStore(const Mesh &mesh, const OutputConfiguration &configuration);
	static void destroyAsynchronizedStore();
	static bool isStoreNode();
	static bool isComputeNode();

	bool isCollected();
	bool isDistributed();
	bool storeStep(const Step &step);

	void updateMesh();
	void updateSolution(const Step &step);

	ResultStore();
	~ResultStore();

	int storeThreads;
	int storeProcesses;
	int computeProcesses;

protected:
	ResultStoreExecutor *_async;
	ResultStoreExecutor *_direct;

	static ResultStore *_asyncStore; // only one instance can be created
	static async::Dispatcher *_dispatcher;
};

}

#endif /* SRC_OUTPUT_RESULT_RESULTSTORE_H_ */
