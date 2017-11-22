
#ifndef SRC_OUTPUT_RESULTS_EXECUTOR_ASYNCEXECUTOR_H_
#define SRC_OUTPUT_RESULTS_EXECUTOR_ASYNCEXECUTOR_H_

#include "directexecutor.h"
#include "async/Module.h"

namespace espreso {

class AsyncBufferManager {

public:
	enum Buffer: int {
		NAME,
		COORDINATES,
		ELEMENTS,
		SOLUTION,

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

	template <typename Buffer, typename...Rest>
	ExecParameters(Buffer buffer, Rest... rest): ExecParameters(rest...) { updatedBuffers |= 1 << buffer; }
};

class AsyncExecutor: public DirectExecutor {

public:
	AsyncExecutor(const OutputConfiguration &configuration): DirectExecutor(configuration) {}

	void execInit(const async::ExecInfo &info, const InitParameters &initParameters);
	void exec(const async::ExecInfo &info, const ExecParameters &parameters);
//	void finalize();
};

class AsyncStore : public ResultStoreExecutor, private async::Module<AsyncExecutor, InitParameters, ExecParameters> {

public:
	AsyncStore(const OutputConfiguration &configuration);
	~AsyncStore();

	void updateMesh();
	void updateSolution();

protected:
	void init();
	void setUp() { setExecutor(_executor); };
//	void tearDown() { _executor.finalize(); }

	AsyncExecutor _executor;
};

}



#endif /* SRC_OUTPUT_RESULTS_EXECUTOR_ASYNCEXECUTOR_H_ */
