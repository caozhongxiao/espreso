
#include "sharedinterfacestore.h"

#include "../../basis/containers/serializededata.h"

using namespace espreso;

SharedInterfaceStore::SharedInterfaceStore()
: nodes(NULL)
{

}

SharedInterfaceStore::~SharedInterfaceStore()
{
	if (nodes == NULL) { delete nodes; }
}



