
#include "fetidatastore.h"

#include "../../basis/containers/serializededata.h"

using namespace espreso;

FETIDataStore::FETIDataStore()
: interfaceNodes(NULL),
  cornerDomains(NULL)
{

}

FETIDataStore::~FETIDataStore()
{
	if (interfaceNodes == NULL) { delete interfaceNodes; }
	if (cornerDomains == NULL) { delete cornerDomains; }
}



