
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

bool ResultStoreExecutor::storeSolution(const Step &step)
{
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
}

bool ResultStoreExecutor::storeStatistics(const Step &step)
{
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
}




