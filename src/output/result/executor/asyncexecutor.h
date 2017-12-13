
#ifndef SRC_OUTPUT_RESULT_EXECUTOR_ASYNCEXECUTOR_H_
#define SRC_OUTPUT_RESULT_EXECUTOR_ASYNCEXECUTOR_H_

#include "async/Module.h"

#include "../../../mesh/mesh.h"
#include "directexecutor.h"

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
	AsyncExecutor(const ECFConfiguration &configuration);

	void execInit(const async::ExecInfo &info, const InitParameters &initParameters);
	void exec(const async::ExecInfo &info, const ExecParameters &parameters);

	const Mesh& mesh() const { return _mesh; }

protected:
	Mesh _mesh;
	const char *_buffer;
};

class AsyncStore: public ResultStoreExecutor, private async::Module<AsyncExecutor, InitParameters, ExecParameters> {

public:
	AsyncStore(const Mesh &mesh, const OutputConfiguration &configuration);
	~AsyncStore();

	bool isCollected() { return true; }
	virtual const Mesh& mesh() const { return _executor.mesh(); }

	void addResultStore(ResultStoreBase *resultStore);
	bool hasStore();

	void updateMesh();
	void updateSolution(const Step &step);

protected:
	void init();
	void setUp() { setExecutor(_executor); };

	void prepareBuffer(AsyncBufferManager::Buffer buffer, size_t size);

	AsyncExecutor _executor;
	char *_buffer;
};

}



#endif /* SRC_OUTPUT_RESULT_EXECUTOR_ASYNCEXECUTOR_H_ */
