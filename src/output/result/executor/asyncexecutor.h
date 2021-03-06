
#ifndef SRC_OUTPUT_RESULT_EXECUTOR_ASYNCEXECUTOR_H_
#define SRC_OUTPUT_RESULT_EXECUTOR_ASYNCEXECUTOR_H_

#include "directexecutor.h"

#ifdef WITHOUT_ASYNC

namespace espreso {
struct AsyncStore: public DirectExecutor {
	AsyncStore(const Mesh &mesh): DirectExecutor(mesh) {}
};
}

#else

//#define USE_MPI
#include "async/Module.h"

#include "mesh/mesh.h"

namespace espreso {

class AsyncBufferManager {

public:
	enum Buffer: int {
		NODES,
		ELEMENTS,

		ELEMENTREGIONS,
		BOUNDARYREGIONS,

		NODEDATA,
		ELEMENTDATA,

		SIZE
	};

	static void buffer(Buffer buffer, int index) { _indices[buffer] = index; }
	static int buffer(Buffer buffer) { return _indices[buffer]; }

protected:
	static std::vector<int> _indices;
};

struct InitParameters
{ };

struct ExecParameters
{
	int updatedBuffers;

	ExecParameters(): updatedBuffers(0) { }

	template <typename...Rest>
	ExecParameters(AsyncBufferManager::Buffer buffer, Rest... rest): ExecParameters(rest...) { updatedBuffers |= 1 << buffer; }
};

class AsyncExecutor: public DirectExecutor {

public:
	AsyncExecutor(ResultStore *store);

	void execInit(const async::ExecInfo &info, const InitParameters &initParameters);
	void exec(const async::ExecInfo &info, const ExecParameters &parameters);

	const Mesh& mesh() const { return _mesh; }

protected:
	Mesh _mesh;
	const char *_buffer;
};

class AsyncStore: public ResultStoreExecutor, private async::Module<AsyncExecutor, InitParameters, ExecParameters> {

public:
	AsyncStore(const Mesh &mesh);
	~AsyncStore();

	virtual const Mesh& mesh() const { return _executor.mesh(); }

	bool isCollected();
	bool isSeparated();

	void addResultStore(ResultStoreBase *resultStore);
	bool hasStore();

	void updateMesh();
	void updateSolution();
	void setParallelStoring(int size, double init, double step);

protected:
	void init();
	void setUp() { setExecutor(_executor); };

	void prepareBuffer(AsyncBufferManager::Buffer buffer, size_t size);

	AsyncExecutor _executor;
	char *_buffer;
};

}

#endif /* WITHOUT_ASYNC */
#endif /* SRC_OUTPUT_RESULT_EXECUTOR_ASYNCEXECUTOR_H_ */
