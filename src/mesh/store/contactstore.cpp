
#include "contactstore.h"

#include "../../basis/containers/serializededata.h"

using namespace espreso;

ContactStore::ContactStore(SurfaceStore *surface)
: eps(0.01), groupsize(1),
  surface(surface),
  xsize(0), ysize(0), zsize(0),
  xbegin(0), xend(0), ybegin(0), yend(0), zbegin(0), zend(0),
  grid(NULL)
{

}

ContactStore::~ContactStore()
{
	if (grid != NULL) { delete grid; }
}


