
#include "elementsregionstore.h"

#include "surfacestore.h"

#include "mesh/elements/element.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"

using namespace espreso;

ElementsRegionStore::ElementsRegionStore(const std::string &name)
: name(name),

  elements(NULL),
  uniqueElements(NULL),

  nodes(NULL),

  uniqueOffset(0),
  uniqueSize(0),
  uniqueTotalSize(0),

  ecounters(static_cast<int>(Element::CODE::SIZE)),

  surface(NULL)
{

}

ElementsRegionStore::~ElementsRegionStore()
{
	if (uniqueElements != NULL && uniqueElements != elements) {
		delete uniqueElements;
	}
	if (elements != NULL) { delete elements; }
	if (nodes != NULL) { delete nodes; }
	if (surface != NULL) { delete surface; }
}

size_t ElementsRegionStore::packedSize() const
{
	if (elements == NULL) {
		return 0;
	}
	return
			utils::packedSize(name) +
			utils::packedSize(uniqueOffset) +
			utils::packedSize(uniqueSize) +
			utils::packedSize(uniqueTotalSize) +
			elements->packedSize() + utils::packedSize(eintervals) +
			nodes->packedSize() +
			utils::packedSize(nintervals) +
			utils::packedSize(ecounters);
}

void ElementsRegionStore::pack(char* &p) const
{
	if (elements == NULL) {
		return;
	}
	utils::pack(name, p);
	utils::pack(uniqueOffset, p);
	utils::pack(uniqueSize, p);
	utils::pack(uniqueTotalSize, p);
	elements->pack(p);
	utils::pack(eintervals, p);
	nodes->pack(p);
	utils::pack(nintervals, p);
	utils::pack(ecounters, p);
}

void ElementsRegionStore::unpack(const char* &p)
{
	utils::unpack(name, p);
	utils::unpack(uniqueOffset, p);
	utils::unpack(uniqueSize, p);
	utils::unpack(uniqueTotalSize, p);
	if (elements == NULL) {
		elements = new serializededata<esint, esint>(1, tarray<esint>(1, 0));
	}
	elements->unpack(p);
	utils::unpack(eintervals, p);
	if (nodes == NULL) {
		nodes = new serializededata<esint, esint>(1, tarray<esint>(1, 0));
	}
	nodes->unpack(p);
	utils::unpack(nintervals, p);
	utils::unpack(ecounters, p);
}

