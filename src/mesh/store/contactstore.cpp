
#include "contactstore.h"

#include "../../basis/containers/serializededata.h"

using namespace espreso;

ContactStore::ContactStore(SurfaceStore *surface)
: eps(0.01), groupsize(1),
  surface(surface),
  elements(NULL),
  closeElements(NULL),
  grid(NULL)
{

}

ContactStore::~ContactStore()
{
	if (grid != NULL) { delete grid; }
	if (elements != NULL) { delete elements; }
	if (closeElements != NULL) { delete closeElements; }
}


