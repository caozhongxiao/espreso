
#include "element.h"

#include "basis/matrices/denseMatrix.h"
#include "basis/containers/serializededata.h"

using namespace espreso;

Element::~Element()
{
	if (N != NULL) { delete N; }
	if (dN != NULL) { delete dN; }
	if (weighFactor != NULL) { delete weighFactor; }

	if (faces != NULL) { delete faces; }
	if (edges != NULL) { delete edges; }
	if (facepointers != NULL) { delete facepointers; }
	if (edgepointers != NULL) { delete edgepointers; }
	if (triangles != NULL) { delete triangles; }
}



