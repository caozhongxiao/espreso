
#include "boundaryregionstore.h"

#include "../../basis/containers/serializededata.h"

using namespace espreso;

BoundaryRegionStore::BoundaryRegionStore(const std::string &name)
: name(name),

  faces(NULL),
  edges(NULL),
  nodes(NULL),

  facepointers(NULL),
  edgepointers(NULL)
{

}

BoundaryRegionStore::~BoundaryRegionStore()
{
	if (faces == NULL) { delete faces; }
	if (edges == NULL) { delete edges; }
	if (nodes == NULL) { delete nodes; }

	if (facepointers == NULL) { delete facepointers; }
	if (edgepointers == NULL) { delete edgepointers; }
}

