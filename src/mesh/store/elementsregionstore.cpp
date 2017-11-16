
#include "elementsregionstore.h"

#include "../../basis/containers/serializededata.h"

using namespace espreso;

ElementsRegionStore::ElementsRegionStore(const std::string &name)
: name(name),

  elements(NULL)
{

}

ElementsRegionStore::~ElementsRegionStore()
{
	if (elements == NULL) { delete elements; }
}



