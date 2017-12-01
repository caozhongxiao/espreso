
#include "elementsregionstore.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/utilities/utils.h"

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

size_t ElementsRegionStore::packedSize() const
{
	return Esutils::packedSize(name) + elements->packedSize();
}

void ElementsRegionStore::pack(char* &p) const
{
	Esutils::pack(name, p);
	elements->pack(p);
}

void ElementsRegionStore::unpack(const char* &p)
{
	Esutils::unpack(name, p);
	if (elements == NULL) {
		elements = new serializededata<eslocal, eslocal>(tarray<eslocal>(1, 0), tarray<eslocal>(1, 0));
	}
	elements->unpack(p);
}

