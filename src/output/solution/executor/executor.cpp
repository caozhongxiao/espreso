
#include "executor.h"

#include "../visualization/vtklegacy.h"
#include "../../../basis/utilities/utils.h"
#include "../../../mesh/mesh.h"

using namespace espreso;

void SolutionStoreExecutor::storePreprocessedData()
{
	std::string root = Esutils::createDirectory({ Logging::outputRoot(), "PREPROCESSED_DATA" });

	VTKLegacy::mesh(root + "/mesh", _configuration, _mesh.nodes, _mesh.elements);
	VTKLegacy::nodesIntervals(root + "/nodeIntervals", _mesh.nodes);
}




