
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

  surface(new SurfaceStore())
{

}

ElementsRegionStore::ElementsRegionStore(const char* &packedData)
: ElementsRegionStore("")
{
	unpackFull(packedData);
}

ElementsRegionStore::~ElementsRegionStore()
{
	if (uniqueElements != NULL && uniqueElements != elements) {
		delete uniqueElements;
	}
	if (elements != NULL) { delete elements; }
	if (nodes != NULL) { delete nodes; }
	delete surface;
}

size_t ElementsRegionStore::packedFullSize() const
{
	size_t packedSize = 0;

	packedSize += utils::packedSize(name);

	packedSize += utils::packedSize(elements);
	packedSize += utils::packedSize(uniqueElements);
	packedSize += utils::packedSize(nodes);

	packedSize += utils::packedSize(eintervals);
	packedSize += utils::packedSize(ueintervals);
	packedSize += utils::packedSize(nintervals);

	packedSize += utils::packedSize(uniqueOffset);
	packedSize += utils::packedSize(uniqueSize);
	packedSize += utils::packedSize(uniqueTotalSize);

	packedSize += utils::packedSize(ecounters);

	return packedSize;
}

void ElementsRegionStore::packFull(char* &p) const
{
	utils::pack(name, p);

	utils::pack(elements, p);
	utils::pack(uniqueElements, p);
	utils::pack(nodes, p);

	utils::pack(eintervals, p);
	utils::pack(ueintervals, p);
	utils::pack(nintervals, p);

	utils::pack(uniqueOffset, p);
	utils::pack(uniqueSize, p);
	utils::pack(uniqueTotalSize, p);

	utils::pack(ecounters, p);
}

void ElementsRegionStore::unpackFull(const char* &p)
{
	utils::unpack(name, p);

	utils::unpack(elements, p);
	utils::unpack(uniqueElements, p);
	utils::unpack(nodes, p);

	utils::unpack(eintervals, p);
	utils::unpack(ueintervals, p);
	utils::unpack(nintervals, p);

	utils::unpack(uniqueOffset, p);
	utils::unpack(uniqueSize, p);
	utils::unpack(uniqueTotalSize, p);

	utils::unpack(ecounters, p);
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

