
#include "executor.h"

#include "../visualization/vtklegacy.h"

#include "../../../basis/utilities/utils.h"

#include "../../../mesh/mesh.h"

using namespace espreso;

void ResultStoreExecutor::storePreprocessedData(const Mesh &mesh)
{
	std::string root = Esutils::createDirectory({ Logging::outputRoot(), "PREPROCESSED_DATA" });

	VTKLegacy::mesh(root + "/mesh", mesh.nodes, mesh.elements);
	VTKLegacy::nodesIntervals(root + "/nodeIntervals", mesh.nodes);
}




