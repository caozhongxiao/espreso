
#include "surfacestore.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"

using namespace espreso;

SurfaceStore::SurfaceStore()
: triangles(NULL),
  elements(NULL),
  nelements(NULL),
  IDs(NULL),
  neighbors(NULL),
  nodes(NULL),
  coordinates(NULL),
  eoffset(0),
  epointers(NULL)
{

}

SurfaceStore::~SurfaceStore()
{
	if (triangles != NULL && triangles != elements) { delete triangles; }
	if (elements != NULL) { delete elements; }
	if (nelements != NULL) { delete nelements; }
	if (IDs != NULL) { delete IDs; }
	if (neighbors != NULL) { delete neighbors; }
	if (nodes != NULL) { delete nodes; }
	if (coordinates != NULL) { delete coordinates; }
	if (epointers != NULL) { delete epointers; }
}


