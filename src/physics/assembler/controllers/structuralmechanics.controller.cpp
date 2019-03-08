
#include "structuralmechanics.controller.h"

using namespace espreso;

StructuralMechanicsController::StructuralMechanicsController(int dimension, StructuralMechanicsController *previous, StructuralMechanicsLoadStepConfiguration &configuration)
: Controller(dimension),
  _configuration(configuration),
  _kcoordinate(dimension), _kdisplacement(dimension),
  _ktemperature(1), _kinitialTemperature(1), _kacceleration(dimension), _kangularVelocity(3), _kthickness(1),
  _ndisplacement(NULL),
  _ethickness(NULL)
{
	if (previous) {
		_ndisplacement = previous->_ndisplacement;
		_ethickness = previous->_ethickness;
	}
}

NodeData* StructuralMechanicsController::solution()
{
	return _ndisplacement;
}

