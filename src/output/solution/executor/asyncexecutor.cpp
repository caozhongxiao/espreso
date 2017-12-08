
#include "asyncexecutor.h"

#include "../../../basis/utilities/utils.h"
#include "../../../config/ecf/ecf.h"

#include "../../../assembler/step.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/elementsregionstore.h"
#include "../../../mesh/store/boundaryregionstore.h"

#include "../visualization/vtklegacy.h"

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

	{ // ELEMENT REGIONS
		size_t esize = 0;
		for (size_t r = 0; r < _mesh.elementsRegions.size(); r++) {
			esize += _mesh.elementsRegions[r]->packedSize();
		}
		prepareBuffer(AsyncBufferManager::ELEMENTREGIONS, esize);

		_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTREGIONS));
		for (size_t r = 0; r < _mesh.elementsRegions.size(); r++) {
			_mesh.elementsRegions[r]->pack(_buffer);
		}
	}

	{ // BOUNDARY REGIONS
		size_t bsize = 0;
		for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
			bsize += _mesh.boundaryRegions[r]->packedSize();
		}
		prepareBuffer(AsyncBufferManager::BOUNDARYREGIONS, bsize);

		_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::BOUNDARYREGIONS));
		for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
			_mesh.boundaryRegions[r]->pack(_buffer);
		}
	}

	call(ExecParameters(AsyncBufferManager::NODES, AsyncBufferManager::ELEMENTS, AsyncBufferManager::ELEMENTREGIONS, AsyncBufferManager::BOUNDARYREGIONS));
}

void AsyncStore::updateSolution(const Step &step)
{
	std::string root = Esutils::createDirectory({ Logging::outputRoot(), "PREPROCESSED_DATA" });

	VTKLegacy::solution(root + "/solution", _configuration, _mesh.nodes, _mesh.elements);

	wait();

	prepareBuffer(AsyncBufferManager::NODEDATA, sizeof(Step) + _mesh.nodes->packedDataSize());
	_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::NODEDATA));
	Esutils::pack(step, _buffer);
	_mesh.nodes->packData(_buffer);

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

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTREGIONS) {
		int bufferid = AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTREGIONS);
		_buffer = static_cast<const char*>(info.buffer(bufferid));
		while (_buffer < info.buffer(bufferid) + info.bufferSize(bufferid)) {
			_mesh.elementsRegions.push_back(new ElementsRegionStore(""));
			_mesh.elementsRegions.back()->unpack(_buffer);
		}
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::BOUNDARYREGIONS) {
		int bufferid = AsyncBufferManager::buffer(AsyncBufferManager::BOUNDARYREGIONS);
		_buffer = static_cast<const char*>(info.buffer(bufferid));
		while (_buffer < info.buffer(bufferid) + info.bufferSize(bufferid)) {
			_mesh.boundaryRegions.push_back(new BoundaryRegionStore("", _mesh._eclasses));
			_mesh.boundaryRegions.back()->unpack(_buffer);
		}
	}

	if (
			(parameters.updatedBuffers & 1 << AsyncBufferManager::NODES) ||
			(parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTS)) {

		updateMesh();
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::NODEDATA) {
		_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::NODEDATA)));
		Step step;
		Esutils::unpack(step, _buffer);
		_mesh.nodes->unpackData(_buffer);
		updateSolution(step);
	}
}

//void AsyncExecutor::finalize()
//{
//	std::cout << "finalize\n";
//}
