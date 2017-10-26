
#include "boundarystore.h"
#include "../../basis/containers/serializededata.h"

using namespace espreso;

BoundaryStore::BoundaryStore()
: elems(NULL),
  faces(NULL),
  edges(NULL),
  nodes(NULL),

  facepointers(),
  edgepointers()
{

}

BoundaryStore::~BoundaryStore()
{
	if (elems == NULL) { delete elems; }
	if (faces == NULL) { delete faces; }
	if (edges == NULL) { delete edges; }
	if (nodes == NULL) { delete nodes; }

	if (facepointers == NULL) { delete facepointers; }
	if (edgepointers == NULL) { delete edgepointers; }
}

