
#include "executor.h"

#include "../../../basis/utilities/utils.h"
#include "../../../config/ecf/output.h"
#include "../../../assembler/step.h"
#include "../../../mesh/mesh.h"

#include "../visualization/distributedvtklegacy.h"
#include "../monitors/monitoring.h"

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

bool ResultStoreExecutor::isDistributed()
{
	for (size_t i = 0; i < _resultStore.size(); i++) {
		if (_resultStore[i]->isDistributed()) {
			return true;
		}
	}
	return false;
}

bool ResultStoreExecutor::storeStep(const Step &step)
{
	return Monitoring::storeStep(_configuration, step) || Visualization::storeStep(_configuration, step);
}



