
#include "mesh/elements/element.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/denseMatrix.h"

namespace espreso {

template<>
void Element::set<Element::CODE::POINT1>()
{
	type = Element::TYPE::POINT;
	code = Element::CODE::POINT1;
	nodes = 1;
	coarseNodes = 1;
	nCommonFace = 1;
	nCommonEdge = 1;
}
}


