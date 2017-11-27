
#include "collectedvisualization.h"

#include "../../../config/ecf/environment.h"

using namespace espreso;

CollectedVisualization::CollectedVisualization()
{
	MPI_Comm_split(environment->MPICommunicator, 0, environment->MPIrank, &_storeCommunicator);
}

CollectedVisualization::~CollectedVisualization()
{
	MPI_Comm_free(&_storeCommunicator);
}


