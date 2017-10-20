
#include "domainstore.h"

#include "../../basis/containers/serializededata.h"

using namespace espreso;

DomainStore::DomainStore()
: size(0),
  elems(NULL),
  nodes(NULL)
{

}

DomainStore::~DomainStore()
{
	if (elems == NULL) { delete elems; }
	if (nodes == NULL) { delete nodes; }
}



