
#include "executor.h"

#include "../visualization/vtklegacy.h"
#include "../../../basis/utilities/utils.h"
#include "../../../config/ecf/output.h"
#include "../../../assembler/step.h"
#include "../../../mesh/mesh.h"

using namespace espreso;

void SolutionStoreExecutor::storePreprocessedData()
{
	std::string root = Esutils::createDirectory({ Logging::outputRoot(), "PREPROCESSED_DATA" });

	VTKLegacy::mesh(root + "/mesh", _configuration, _mesh.nodes, _mesh.elements);
	VTKLegacy::nodesIntervals(root + "/nodeIntervals", _mesh.nodes);
}

bool SolutionStoreExecutor::storeSolution(const Step &step)
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




