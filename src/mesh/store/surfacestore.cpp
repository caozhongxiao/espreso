
#include "surfacestore.h"

#include "../../basis/containers/point.h"
#include "../../basis/containers/serializededata.h"

using namespace espreso;

SurfaceStore::SurfaceStore()
: triangles(NULL), elements(NULL), coordinates(NULL), epointers(NULL)
{

}

SurfaceStore::~SurfaceStore()
{
	if (triangles != NULL && triangles != elements) { delete triangles; }
	if (elements != NULL) { delete elements; }
	if (coordinates != NULL) { delete coordinates; }
	if (epointers != NULL) { delete epointers; }
}


