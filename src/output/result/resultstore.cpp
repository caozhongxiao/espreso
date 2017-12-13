
#include "../result/resultstore.h"

#include "../../config/ecf/environment.h"
#include "../../config/ecf/output.h"
#include "../../basis/logging/logging.h"
#include "executor/asyncexecutor.h"
#include "executor/directexecutor.h"
#include "async/Dispatcher.h"

#include "monitors/monitoring.h"
#include "visualization/collectedensight.h"


using namespace espreso;

ResultStore* ResultStore::_asyncStore = NULL;
async::Dispatcher* ResultStore::_dispatcher = NULL;

ResultStore* ResultStore::createAsynchronizedStore(const Mesh &mesh, const OutputConfiguration &configuration)
{
	if (_asyncStore != NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: try to create multiple ASYNC instances.";
	}
	_asyncStore = new ResultStore();

	_asyncStore->_direct = new DirectExecutor(mesh, configuration);
	switch (configuration.mode) {
	case OutputConfiguration::MODE::SYNC:
		break;
	case OutputConfiguration::MODE::THREAD:
		_dispatcher = new async::Dispatcher();
		_asyncStore->_async = new AsyncStore(mesh, configuration);
		async::Config::setMode(async::THREAD);
		break;
	case OutputConfiguration::MODE::MPI:
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: not implemented OUTPUT::MODE==MPI.";
		_dispatcher = new async::Dispatcher();
		_asyncStore->_async = new AsyncStore(mesh, configuration);
		if (environment->MPIsize == 1) {
			ESINFO(GLOBAL_ERROR) << "Invalid number of MPI processes. OUTPUT::MODE==MPI required at least two MPI processes.";
		}
		if (configuration.output_node_group_size == 0) {
			ESINFO(GLOBAL_ERROR) << "OUTPUT::OUTPUT_NODE_GROUP_SIZE cannot be 0.";
		}
		async::Config::setMode(async::MPI);
		async::Config::setGroupSize(configuration.output_node_group_size);
		_dispatcher->setGroupSize(configuration.output_node_group_size);
		break;
	}

	// TODO: optimize
	_asyncStore->_async->addResultStore(new CollectedEnSightWithDecomposition("solution", _asyncStore->_async->mesh()));
	if (configuration.monitoring.size()) {
		_asyncStore->_async->addResultStore(new Monitoring(_asyncStore->_async->mesh(), configuration, true));
		// _asyncStore->_direct->addResultStore(new Monitoring(_asyncStore->_direct->mesh(), configuration, false));
	}

	if (_asyncStore->_direct->hasStore() == 0) {
		delete _asyncStore->_direct;
		_asyncStore->_direct = NULL;
	}

	if (_asyncStore->_async->hasStore() == 0) {
		delete _asyncStore->_dispatcher;
		delete _asyncStore->_async;
		_asyncStore->_dispatcher = NULL;
		_asyncStore->_async = NULL;
	}

	if (_dispatcher != NULL) {
		_dispatcher->init();
		_dispatcher->dispatch();

		if (configuration.mode == OutputConfiguration::MODE::MPI) {
			int computeSize;
			MPI_Comm_size(_dispatcher->commWorld(), &computeSize);
			ESINFO(OVERVIEW) << "Reserve " << environment->MPIsize - computeSize << " MPI process(es) for storing solution.";
		}

		environment->MPICommunicator = _dispatcher->commWorld();
		MPI_Comm_rank(environment->MPICommunicator, &environment->MPIrank);
		MPI_Comm_size(environment->MPICommunicator, &environment->MPIsize);
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

bool ResultStore::storeStep(const Step &step)
{
	bool store = false;
	if (_async) store |= _async->storeStep(step);
	if (_direct) store |= _direct->storeStep(step);
	return store;
}

void ResultStore::updateMesh()
{
	if (_async) _async->updateMesh();
	if (_direct) _direct->updateMesh();
}

void ResultStore::updateSolution(const Step &step)
{
	if (_async) _async->updateSolution(step);
	if (_direct) _direct->updateSolution(step);
}

ResultStore::ResultStore(): _async(NULL), _direct(NULL)
{

}

ResultStore::~ResultStore()
{
	if (_async) delete _async;
	if (_direct) delete _direct;
}


