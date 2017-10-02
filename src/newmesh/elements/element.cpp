
#include "element.h"

#include "../../basis/containers/serializededata.h"
#include "elementstore.h"

using namespace espreso;

edata<eslocal> NewElement::nodes()
{
	return *(_store->nodes->begin() + _index);
}

edata<const eslocal> NewElement::nodes() const
{
	return *(_store->nodes->cbegin() + _index);
}



