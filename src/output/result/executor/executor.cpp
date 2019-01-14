
#include "executor.h"

#include "basis/utilities/utils.h"
#include "config/ecf/output.h"
#include "mesh/mesh.h"

#include "output/result/monitors/monitoring.h"
#include "output/result/visualization/separated/vtklegacy.h"

using namespace espreso;

ResultStoreExecutor::~ResultStoreExecutor()
{
	for (size_t i = 0; i < _resultStore.size(); i++) {
		delete _resultStore[i];
	}
}

void ResultStoreExecutor::addResultStore(ResultStoreBase *resultStore)
{
	_resultStore.push_back(resultStore);
}

bool ResultStoreExecutor::isCollected()
{
	for (size_t i = 0; i < _resultStore.size(); i++) {
		if (_resultStore[i]->isCollected()) {
			return true;
		}
	}
	return false;
}

bool ResultStoreExecutor::isSeparated()
{
	for (size_t i = 0; i < _resultStore.size(); i++) {
		if (_resultStore[i]->isSeparated()) {
			return true;
		}
	}
	return false;
}

bool ResultStoreExecutor::storeStep()
{
	return Monitoring::storeStep(_configuration) || Visualization::storeStep(_configuration);
}




