
#include "elementwrapper.h"

#include "../../basis/containers/serializededata.h"
#include "elementstore.h"

using namespace espreso;

edata<eslocal> ElementWrapper::nodes()
{
	return *(_store->nodes->begin() + _index);
}

edata<const eslocal> ElementWrapper::nodes() const
{
	return *(_store->nodes->cbegin() + _index);
}



