
#include "resultstore.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/output.h"
#include "basis/utilities/sysutils.h"

#include "executor/asyncexecutor.h"
#include "executor/directexecutor.h"

#include "monitors/monitoring.h"
#include "visualization/collected/ensight.h"
#include "visualization/collected/stl.h"
#include "visualization/separated/insitu.h"
#include "visualization/separated/vtklegacy.h"

#include <string>

#ifdef WITHOUT_ASYNC
struct espreso::Dispatcher {
	void init() {}
	void dispatch() {}
	bool isExecutor() { return false; }
	MPI_Comm commWorld() { return MPI_COMM_WORLD; }
	void setThreadMode() {}
};
#else
#include "async/Dispatcher.h"
struct espreso::Dispatcher: public async::Dispatcher {
void setThreadMode() { async::Config::setMode(async::THREAD); }
};
#endif /* WITHOUT_ASYNC */


using namespace espreso;

ResultStore* ResultStore::_asyncStore = NULL;
Dispatcher* ResultStore::_dispatcher = NULL;

ResultStoreBase::ResultStoreBase(const Mesh &mesh): _mesh(mesh), _directory("PREPOSTDATA/")
{
}

void ResultStoreBase::createOutputDirectory()
{
	utils::createDirectory({ std::string(eslog::path()), _directory });
}

ResultStore* ResultStore::createAsynchronizedStore(const Mesh &mesh)
{
	if (_asyncStore != NULL) {
		eslog::globalerror("ESPRESO internal error: try to create multiple ASYNC instances.\n");
	}


	_asyncStore = new ResultStore();
	_asyncStore->storeThreads = 0;
	_asyncStore->storeProcesses = 0;
	_asyncStore->computeProcesses = info::mpi::size;

	_asyncStore->_direct = new DirectExecutor(mesh);
	ResultStoreExecutor *executor = _asyncStore->_direct;
	switch (info::ecf->output.mode) {
	case OutputConfiguration::MODE::SYNC:
		break;
	case OutputConfiguration::MODE::THREAD:
		_dispatcher = new Dispatcher();
		_asyncStore->_async = new AsyncStore(mesh);
		_dispatcher->setThreadMode();
		_asyncStore->storeThreads = 1;
		executor = _asyncStore->_async;
		break;
//	case OutputConfiguration::MODE::MPI:
//		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: not implemented OUTPUT::MODE==MPI.";
//		_dispatcher = new async::Dispatcher();
//		_asyncStore->_async = new AsyncStore(mesh, configuration);
//		if (info::mpi::MPIsize == 1) {
//			ESINFO(GLOBAL_ERROR) << "Invalid number of MPI processes. OUTPUT::MODE==MPI required at least two MPI processes.";
//		}
//		if (configuration.output_node_group_size == 0) {
//			ESINFO(GLOBAL_ERROR) << "OUTPUT::OUTPUT_NODE_GROUP_SIZE cannot be 0.";
//		}
//		async::Config::setMode(async::MPI);
//		async::Config::setGroupSize(configuration.output_node_group_size);
//		_dispatcher->setGroupSize(configuration.output_node_group_size);
//		break;
	}

	// TODO: optimize
	if (info::ecf->output.results_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
		switch (info::ecf->output.format) {
		case OutputConfiguration::FORMAT::ENSIGHT:
			executor->addResultStore(new EnSightWithDecomposition(std::string(eslog::name()), executor->mesh()));
			break;
		case OutputConfiguration::FORMAT::STL_SURFACE:
			executor->addResultStore(new STL(std::string(eslog::name()), mesh));
			break;
		default:
			eslog::globalerror("ESPRESO internal error: implement the selected output format.\n");
		}

	}
	if (info::ecf->output.monitors_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER && info::ecf->output.monitoring.size()) {
		_asyncStore->_direct->addResultStore(new Monitoring(mesh));
	}
	if (info::ecf->output.catalyst) {
		_asyncStore->_direct->addResultStore(new InSitu(mesh));
	}
	if (info::ecf->output.debug) {
		_asyncStore->_direct->addResultStore(new VTKLegacyDebugInfo(mesh));
	}

	if (!_asyncStore->_direct->hasStore()) {
		delete _asyncStore->_direct;
		_asyncStore->_direct = NULL;
	}

	if (_asyncStore->_async != NULL && !_asyncStore->_async->hasStore()) {
		delete _asyncStore->_dispatcher;
		delete _asyncStore->_async;
		_asyncStore->_dispatcher = NULL;
		_asyncStore->_async = NULL;
	}

	if (_dispatcher != NULL) {
		_dispatcher->init();
		_dispatcher->dispatch();

//		if (configuration.mode == OutputConfiguration::MODE::MPI) {
//			int computeSize;
//			MPI_Comm_size(_dispatcher->commWorld(), &computeSize);
//			_asyncStore->computeProcesses = computeSize;
//			_asyncStore->storeProcesses = info::mpi::MPIsize - computeSize;
//		}

		info::mpi::comm = _dispatcher->commWorld();
		MPI_Comm_rank(info::mpi::comm, &info::mpi::rank);
		MPI_Comm_size(info::mpi::comm, &info::mpi::size);
	}

	return _asyncStore;
}

void ResultStore::destroyAsynchronizedStore()
{
	if (_asyncStore) {
		delete _asyncStore;
	}
	if (_dispatcher) {
		delete _dispatcher;
	}

	_asyncStore = NULL;
	_dispatcher = NULL;
}

bool ResultStore::isStoreNode()
{
	return _dispatcher != NULL && _dispatcher->isExecutor();
}

bool ResultStore::isComputeNode()
{
	return !isStoreNode();
}

bool ResultStore::isCollected()
{
	bool store = false;
	if (_async) store |= _async->isCollected();
	if (_direct) store |= _direct->isCollected();
	return store;
}

bool ResultStore::isSeparated()
{
	bool store = false;
	if (_async) store |= _async->isSeparated();
	if (_direct) store |= _direct->isSeparated();
	return store;
}

bool ResultStore::storeStep()
{
	bool store = false;
	if (_async) store |= _async->storeStep();
	if (_direct) store |= _direct->storeStep();
	return store;
}

void ResultStore::updateMesh()
{
	if (_async && _async->hasStore()) _async->updateMesh();
	if (_direct && _direct->hasStore()) _direct->updateMesh();
}

void ResultStore::updateSolution()
{
	if (_async && storeStep()) _async->updateSolution();
	if (_direct && storeStep()) _direct->updateSolution();
}

ResultStore::ResultStore(): _async(NULL), _direct(NULL)
{

}

ResultStore::~ResultStore()
{
	if (_async) delete _async;
	if (_direct) delete _direct;
}


