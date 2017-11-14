
#include "domainstore.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/logging/logging.h"

using namespace espreso;

DomainStore::DomainStore()
: offset(0),
  size(0),
  elems(NULL),
  nodes(NULL)
{

}

DomainStore::~DomainStore()
{
	if (elems == NULL) { delete elems; }
	if (nodes == NULL) { delete nodes; }
}

std::vector<esglobal> DomainStore::gatherElementDistribution()
{
	esglobal offset = domainElementBoundaries.back();
	Communication::exscan(offset);

	std::vector<esglobal> distribution(domainElementBoundaries.begin(), domainElementBoundaries.end());
	for (size_t i = 0; i < distribution.size(); i++) {
		distribution[i] += offset;
	}
	std::vector<esglobal> result;

	if (!Communication::gatherUnknownSize(distribution, result)) {
		ESINFO(ERROR) << "ESPRESO internal error: gather domain distribution.";
	}
	if (!Communication::broadcastUnknownSize(result)) {
		ESINFO(ERROR) << "ESPRESO internal error: broadcast domain distribution.";
	}

	Esutils::removeDuplicity(result);
	return result;
}

std::vector<esglobal> DomainStore::gatherProcsDistribution()
{
	std::vector<esglobal> result(environment->MPIsize + 1);

	MPI_Allgather(&offset, sizeof(eslocal), MPI_BYTE, result.data(), sizeof(eslocal), MPI_BYTE, environment->MPICommunicator);
	result.back() = offset + size;
	MPI_Bcast(&result.back(), sizeof(eslocal), MPI_BYTE, environment->MPIsize - 1, environment->MPICommunicator);

	return result;
}

