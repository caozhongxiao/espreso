
#include "element.h"

#include "../../basis/containers/serializededata.h"
#include "elementstore.h"


using namespace espreso;

edata<eslocal> NewElement::indices()
{
	return *(_store->indices->begin() + _index);
}

edata<const eslocal> NewElement::indices() const
{
	return *(_store->indices->cbegin() + _index);
}



