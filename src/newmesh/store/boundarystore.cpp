
#include "boundarystore.h"
#include "../../basis/containers/serializededata.h"

using namespace espreso;

BoundaryStore::BoundaryStore()
: clusterfaces(NULL),
  localfaces(NULL),
  edges(NULL),
  nodes(NULL),

  facepointers(),
  edgepointers()
{

}

BoundaryStore::~BoundaryStore()
{
	if (clusterfaces == NULL) { delete clusterfaces; }
	if (localfaces == NULL) { delete localfaces; }
	if (edges == NULL) { delete edges; }
	if (nodes == NULL) { delete nodes; }

	if (facepointers == NULL) { delete facepointers; }
	if (edgepointers == NULL) { delete edgepointers; }
}



