
#ifndef SRC_OUTPUT_RESULT_RESULTSTORE_H_
#define SRC_OUTPUT_RESULT_RESULTSTORE_H_

namespace async { class Dispatcher; }

namespace espreso {

struct Step;
class Mesh;
class OutputConfiguration;
class ResultStoreExecutor;

class ResultStoreBase {

public:
	virtual bool isCollected() =0;

	virtual void updateMesh() =0;
	virtual void updateSolution(const Step &step) =0;

	virtual ~ResultStoreBase() {};

protected:
	ResultStoreBase(const Mesh &mesh): _mesh(mesh) {}

	const Mesh &_mesh;
};

class ResultStore {

public:
	static ResultStore* createAsynchronizedStore(const Mesh &mesh, const OutputConfiguration &configuration);
	static void destroyAsynchronizedStore();
	static bool isStoreNode();
	static bool isComputeNode();

	bool isCollected();
	bool storeStep(const Step &step);

	void updateMesh();
	void updateSolution(const Step &step);

	ResultStore();
	~ResultStore();

protected:
	ResultStoreExecutor *_async;
	ResultStoreExecutor *_direct;

	static ResultStore *_asyncStore; // only one instance can be created
	static async::Dispatcher *_dispatcher;
};

}

#endif /* SRC_OUTPUT_RESULT_RESULTSTORE_H_ */
