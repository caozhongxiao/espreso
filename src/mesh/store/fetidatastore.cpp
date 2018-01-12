
#include "fetidatastore.h"

#include "../../basis/containers/serializededata.h"

using namespace espreso;

FETIDataStore::FETIDataStore()
: interfaceNodes(NULL),
  cornerDomains(NULL)
{
	sFixPointsDistribution = {0, 0};
	iFixPointsDistribution = {0, 0};
}

FETIDataStore::~FETIDataStore()
{
	if (interfaceNodes == NULL) { delete interfaceNodes; }
	if (cornerDomains == NULL) { delete cornerDomains; }
}



