
#ifndef SRC_OUTPUT_RESULT_RESULTSTORE_H_
#define SRC_OUTPUT_RESULT_RESULTSTORE_H_

#include <string>

namespace espreso {

class Mesh;
class ResultStoreExecutor;
class Dispatcher;

class ResultStoreBase {

public:
	virtual bool isCollected() =0;
	virtual bool isSeparated() =0;

	virtual void updateMesh() =0;
	virtual void updateSolution() =0;

	virtual const Mesh& mesh() const { return _mesh; }

	virtual ~ResultStoreBase() {};

protected:
	ResultStoreBase(const Mesh &mesh);

	void createOutputDirectory();

	const Mesh &_mesh;
	std::string _directory;
};

class ResultStore {

public:
	static ResultStore* createAsynchronizedStore(const Mesh &mesh);
	static void destroyAsynchronizedStore();
	static bool isStoreNode();
	static bool isComputeNode();

	bool storeStep();

	void updateMesh();
	void updateSolution();

	ResultStore();
	~ResultStore();

protected:
	ResultStoreExecutor *_async;
	ResultStoreExecutor *_direct;

	static ResultStore *_asyncStore; // only one instance can be created
	static Dispatcher *_dispatcher;
};

}

#endif /* SRC_OUTPUT_RESULT_RESULTSTORE_H_ */
