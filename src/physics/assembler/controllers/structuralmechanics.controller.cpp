
#include "structuralmechanics.controller.h"

#include "mesh/store/nodestore.h"

using namespace espreso;

NodeData* StructuralMechanicsController::solution()
{
	return _displacement;
}

