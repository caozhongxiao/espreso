
#include "contactstore.h"

#include "basis/containers/serializededata.h"

using namespace espreso;

ContactStore::ContactStore(SurfaceStore *surface)
: eps(0.01), groupsize(1),
  surface(surface),
  elements(NULL),
  enormals(NULL),
  closeElements(NULL),
  grid(NULL)
{

}

ContactStore::~ContactStore()
{
	if (grid != NULL) { delete grid; }
	if (elements != NULL) { delete elements; }
	if (enormals != NULL) { delete enormals; }
	if (closeElements != NULL) { delete closeElements; }

	for (size_t i = 0; i < gnelements.size(); i++) {
		delete gnelements[i];
	}
	for (size_t i = 0; i < gngrid.size(); i++) {
		delete gngrid[i];
	}
	for (size_t i = 0; i < gncloseElements.size(); i++) {
		delete gncloseElements[i];
	}
}

size_t ContactStore::packedFullSize() const
{
	size_t packedSize = 0;
	return packedSize;
}

void ContactStore::packFull(char* &p) const
{

}

void ContactStore::unpackFull(const char* &p)
{

}

