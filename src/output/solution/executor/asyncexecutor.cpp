
#include "asyncexecutor.h"

#include "../../../config/ecf/ecf.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"

#include <iostream>

using namespace espreso;

std::vector<int> AsyncBufferManager::_indices(Buffer::SIZE, -1);

void AsyncStore::prepareBuffer(AsyncBufferManager::Buffer buffer, size_t size)
{
	int index = AsyncBufferManager::buffer(buffer);
	if (index == -1) {
		AsyncBufferManager::buffer(buffer, index = addBuffer(NULL, size));
	}
	if (bufferSize(index) != size) {
		removeBuffer(index);
		AsyncBufferManager::buffer(buffer, index = addBuffer(NULL, size));
	}
}

void AsyncStore::updateMesh()
{
	wait();

	prepareBuffer(AsyncBufferManager::NODES, _mesh.nodes->packedSize());
	_mesh.nodes->pack(_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::NODES)));

	prepareBuffer(AsyncBufferManager::ELEMENTS, _mesh.elements->packedSize());
	_mesh.elements->pack(_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTS)));

	call(ExecParameters(AsyncBufferManager::NODES, AsyncBufferManager::ELEMENTS));
}

void AsyncStore::updateSolution()
{
	wait();

	prepareBuffer(AsyncBufferManager::NODEDATA, _mesh.nodes->packedDataSize());
	_mesh.nodes->packData(_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::NODEDATA)));

	call(ExecParameters(AsyncBufferManager::NODEDATA));
}

AsyncStore::AsyncStore(const Mesh &mesh, const OutputConfiguration &configuration)
: SolutionStoreExecutor(mesh, configuration), _executor(mesh.configuration), _buffer(NULL)
{
	async::Module<AsyncExecutor, InitParameters, ExecParameters>::init();
	callInit(InitParameters());
}

AsyncStore::~AsyncStore()
{
	wait();
	async::Module<AsyncExecutor, InitParameters, ExecParameters>::finalize();
}

AsyncExecutor::AsyncExecutor(const ECFConfiguration &configuration)
: _mesh(configuration), DirectExecutor(_mesh, configuration.output), _buffer(NULL)
{

}

void AsyncExecutor::execInit(const async::ExecInfo &info, const InitParameters &initParameters)
{

}

void AsyncExecutor::exec(const async::ExecInfo &info, const ExecParameters &parameters)
{
	if (parameters.updatedBuffers & 1 << AsyncBufferManager::NODES) {
		_mesh.nodes->unpack(_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::NODES))));
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTS) {
		_mesh.elements->unpack(_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTS))));
	}

	if (
			(parameters.updatedBuffers & 1 << AsyncBufferManager::NODES) ||
			(parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTS)) {

		updateMesh();
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::NODEDATA) {
		_mesh.nodes->unpackData(_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::NODEDATA))));
		updateSolution();
	}
}

//void AsyncExecutor::finalize()
//{
//	std::cout << "finalize\n";
//}
