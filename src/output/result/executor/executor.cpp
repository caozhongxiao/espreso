
#include "executor.h"

#include "../../../basis/utilities/utils.h"
#include "../../../config/ecf/output.h"
#include "../../../assembler/step.h"
#include "../../../mesh/mesh.h"
#include "../visualization/vtklegacy.h"

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
	auto storeVisualization = [&] () {
		switch (_configuration.results_store_frequency) {
		case OutputConfiguration::STORE_FREQUENCY::NEVER:
			return false;
		case OutputConfiguration::STORE_FREQUENCY::EVERY_TIMESTEP:
			return true;
		case OutputConfiguration::STORE_FREQUENCY::EVERY_NTH_TIMESTEP:
			return step.substep % _configuration.results_nth_stepping == 0;
		case OutputConfiguration::STORE_FREQUENCY::LAST_TIMESTEP:
			return step.isLast();
		default:
			return false;
		}
	};

	auto storeMonitoring = [&] () {
		switch (_configuration.monitors_store_frequency) {
		case OutputConfiguration::STORE_FREQUENCY::NEVER:
			return false;
		case OutputConfiguration::STORE_FREQUENCY::EVERY_TIMESTEP:
			return true;
		case OutputConfiguration::STORE_FREQUENCY::EVERY_NTH_TIMESTEP:
			return step.substep % _configuration.monitors_nth_stepping == 0;
		case OutputConfiguration::STORE_FREQUENCY::LAST_TIMESTEP:
			return step.isLast();
		default:
			return false;
		}
	};

	return storeMonitoring() || storeVisualization();
}




