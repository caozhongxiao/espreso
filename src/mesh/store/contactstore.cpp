
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

	for (size_t i = 0; i < nelements.size(); i++) {
		delete nelements[i];
	}
	for (size_t i = 0; i < ngrid.size(); i++) {
		delete ngrid[i];
	}
	for (size_t i = 0; i < ncloseElements.size(); i++) {
		delete ncloseElements[i];
	}
}


