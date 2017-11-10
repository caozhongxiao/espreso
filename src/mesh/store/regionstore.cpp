
#include "regionstore.h"
#include "../../basis/containers/serializededata.h"

using namespace espreso;

RegionStore::RegionStore(const std::string &name, TFlags::ELEVEL etype)
: name(name),
  etype(etype),

  elems(NULL),
  faces(NULL),
  edges(NULL),
  nodes(NULL),

  facepointers(),
  edgepointers()
{

}

RegionStore::~RegionStore()
{
	if (elems == NULL) { delete elems; }
	if (faces == NULL) { delete faces; }
	if (edges == NULL) { delete edges; }
	if (nodes == NULL) { delete nodes; }

	if (facepointers == NULL) { delete facepointers; }
	if (edgepointers == NULL) { delete edgepointers; }
}
