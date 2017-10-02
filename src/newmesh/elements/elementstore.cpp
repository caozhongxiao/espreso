
#include "elementstore.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/point/point.h"

using namespace espreso;

ElementStore::ElementStore()
: IDs(NULL),

  elems(NULL),
  faces(NULL),
  edges(NULL),
  nodes(NULL),

  coordinates(NULL),

  globalDual(NULL),
  localDual(NULL)
{

}

ElementStore::~ElementStore()
{
	if (IDs == NULL) { delete IDs; }

	if (elems == NULL) { delete elems; }
	if (faces == NULL) { delete faces; }
	if (edges == NULL) { delete edges; }
	if (nodes == NULL) { delete nodes; }

	if (coordinates == NULL) { delete coordinates; }

	if (globalDual == NULL) { delete globalDual; }
	if (localDual == NULL) { delete localDual; }
}
