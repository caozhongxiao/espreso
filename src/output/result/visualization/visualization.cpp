
#include "visualization.h"

#include "../../../config/ecf/output.h"
#include "../../../physics/step.h"

using namespace espreso;

Visualization::Visualization(const Mesh &mesh, const OutputConfiguration &configuration)
: ResultStoreBase(mesh), _configuration(configuration)
{
	if (configuration.results_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
		createOutputDirectory();
	}
}

bool Visualization::storeStep(const OutputConfiguration &configuration, const Step &step)
{
	switch (configuration.results_store_frequency) {
	case OutputConfiguration::STORE_FREQUENCY::NEVER:
		return false;
	case OutputConfiguration::STORE_FREQUENCY::EVERY_TIMESTEP:
		return true;
	case OutputConfiguration::STORE_FREQUENCY::EVERY_NTH_TIMESTEP:
		return step.substep % configuration.results_nth_stepping == 0;
	case OutputConfiguration::STORE_FREQUENCY::LAST_TIMESTEP:
		return step.isLast();
	default:
		return false;
	}
}


