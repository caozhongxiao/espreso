
#include "elementstore.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/point/point.h"

using namespace espreso;

ElementStore::ElementStore()
: size(0),
  distribution({0, 0}),
  IDs(NULL),

  elems(NULL),
  faces(NULL),
  edges(NULL),
  nodes(NULL),

  coordinates(NULL),
  body(NULL),
  material(NULL),
  epointers(NULL),
  ranks(NULL),

  dual(NULL),
  decomposedDual(NULL)
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
	if (body == NULL) { delete body; }
	if (material == NULL) { delete material; }
	if (epointers == NULL) { delete epointers; }
	if (ranks == NULL) { delete ranks; }

	if (dual == NULL) { delete dual; }
	if (decomposedDual == NULL) { delete decomposedDual; }
}
