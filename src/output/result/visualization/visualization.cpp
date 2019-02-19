
#include "esinfo/timeinfo.h"
#include "visualization.h"

#include "esinfo/ecfinfo.h"
#include "config/ecf/output.h"

using namespace espreso;

Visualization::Visualization(const Mesh &mesh)
: ResultStoreBase(mesh)
{
	if (info::ecf->output.results_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
		createOutputDirectory();
	}
}

bool Visualization::storeStep()
{
	switch (info::ecf->output.results_store_frequency) {
	case OutputConfiguration::STORE_FREQUENCY::NEVER:
		return false;
	case OutputConfiguration::STORE_FREQUENCY::EVERY_TIMESTEP:
		return true;
	case OutputConfiguration::STORE_FREQUENCY::EVERY_NTH_TIMESTEP:
		return time::substep % info::ecf->output.results_nth_stepping == 0;
	case OutputConfiguration::STORE_FREQUENCY::LAST_TIMESTEP:
		return time::isLast();
	default:
		return false;
	}
}


