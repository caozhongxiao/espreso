
#include "visualization.h"

#include "../../../config/ecf/output.h"
#include "../../../globals/time.h"

using namespace espreso;

Visualization::Visualization(const Mesh &mesh, const OutputConfiguration &configuration)
: ResultStoreBase(mesh), _configuration(configuration)
{
	if (configuration.results_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
		createOutputDirectory();
	}
}

bool Visualization::storeStep(const OutputConfiguration &configuration)
{
	switch (configuration.results_store_frequency) {
	case OutputConfiguration::STORE_FREQUENCY::NEVER:
		return false;
	case OutputConfiguration::STORE_FREQUENCY::EVERY_TIMESTEP:
		return true;
	case OutputConfiguration::STORE_FREQUENCY::EVERY_NTH_TIMESTEP:
		return time::substep % configuration.results_nth_stepping == 0;
	case OutputConfiguration::STORE_FREQUENCY::LAST_TIMESTEP:
		return time::isLast();
	default:
		return false;
	}
}


