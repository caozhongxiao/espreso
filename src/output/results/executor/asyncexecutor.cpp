
#include "asyncexecutor.h"

#include <iostream>

using namespace espreso;

std::vector<int> AsyncBufferManager::_indices(Buffer::SIZE);

void AsyncStore::updateMesh()
{
	wait();

	double x = 0.55;
	memcpy(managedBuffer<char*>(AsyncBufferManager::COORDINATES), &x, sizeof(double));
	sendBuffer(AsyncBufferManager::buffer(AsyncBufferManager::COORDINATES));

	int a = 3542;
	memcpy(managedBuffer<char*>(AsyncBufferManager::ELEMENTS), &a, sizeof(int));
	sendBuffer(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTS));

	call(ExecParameters(AsyncBufferManager::COORDINATES, AsyncBufferManager::ELEMENTS));
}

void AsyncStore::updateSolution()
{
	wait();

	int a = 100;
	memcpy(managedBuffer<char*>(AsyncBufferManager::SOLUTION), &a, sizeof(int));
	sendBuffer(AsyncBufferManager::buffer(AsyncBufferManager::SOLUTION));

	call(ExecParameters(AsyncBufferManager::SOLUTION));
}

AsyncStore::AsyncStore(const OutputConfiguration &configuration): ResultStoreExecutor(configuration), _executor(configuration)
{
	async::Module<AsyncExecutor, InitParameters, ExecParameters>::init();

	AsyncBufferManager::buffer(AsyncBufferManager::NAME, addBuffer(NULL, 1024));
	AsyncBufferManager::buffer(AsyncBufferManager::COORDINATES, addBuffer(NULL, sizeof(double)));
	AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTS, addBuffer(NULL, sizeof(int)));
	AsyncBufferManager::buffer(AsyncBufferManager::SOLUTION, addBuffer(NULL, sizeof(int)));

	callInit(InitParameters());
}

AsyncStore::~AsyncStore()
{
	wait();
	async::Module<AsyncExecutor, InitParameters, ExecParameters>::finalize();
}

void AsyncExecutor::execInit(const async::ExecInfo &info, const InitParameters &initParameters)
{

}

void AsyncExecutor::exec(const async::ExecInfo &info, const ExecParameters &parameters)
{
	if (parameters.updatedBuffers & 1 << AsyncBufferManager::COORDINATES) {
		std::cout << "COORDINATES: " << *static_cast<const double*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::COORDINATES))) << "\n";
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTS) {
		std::cout << "ELEMENTS: " << *static_cast<const int*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTS))) << "\n";
	}

	if (
			(parameters.updatedBuffers & 1 << AsyncBufferManager::COORDINATES) ||
			(parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTS)) {

		updateMesh();
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::SOLUTION) {
		std::cout << "SOLUTION: " << *static_cast<const int*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::SOLUTION))) << "\n";
		updateSolution();
	}
}

//void AsyncExecutor::finalize()
//{
//	std::cout << "finalize\n";
//}
