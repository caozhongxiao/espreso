
#include "elementstore.h"
#include "../../basis/containers/serializededata.h"
#include "../../basis/point/point.h"
#include "../../basis/logging/logging.h"

#include <algorithm>
#include <numeric>

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

void ElementStore::sort()
{
	if (IDs == NULL) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: sort element store that has no IDs.";
	}
	std::vector<eslocal> permutation(IDs->datatarray().size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return IDs->datatarray().data()[i] < IDs->datatarray().data()[j]; });

	permute(permutation);
}

void ElementStore::permute(const std::vector<eslocal> &permutation)
{
	if (IDs != NULL) { IDs->permute(permutation); }

	if (elems != NULL) { elems->permute(permutation); }
	if (faces != NULL) { faces->permute(permutation); }
	if (edges != NULL) { edges->permute(permutation); }
	if (nodes != NULL) { nodes->permute(permutation); }

	if (coordinates != NULL) { coordinates->permute(permutation); }
	if (body != NULL) { body->permute(permutation); }
	if (material != NULL) { material->permute(permutation); }
	if (epointers != NULL) { epointers->permute(permutation); }
	if (ranks != NULL) { ranks->permute(permutation); }

	if (dual != NULL) { dual->permute(permutation); }
	if (decomposedDual != NULL) { decomposedDual->permute(permutation); }
}
