
#include "surfacestore.h"

#include "../../basis/containers/serializededata.h"

using namespace espreso;

SurfaceStore::SurfaceStore()
: triangles(NULL), coordinates(NULL)
{

}

SurfaceStore::~SurfaceStore()
{
	if (triangles == NULL) { delete triangles; }
	if (coordinates == NULL) { delete coordinates; }
}


