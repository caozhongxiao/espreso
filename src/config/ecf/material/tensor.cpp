
#include "tensor.h"

#include "../../../basis/logging/logging.h"
#include <iostream>

using namespace espreso;

TensorConfiguration::TensorConfiguration(size_t size)
: size(size), values(size * size)
{

}

ECFExpression& TensorConfiguration::value(size_t row, size_t column)
{
	if (row < size && column < size) {
		return values[row * size + column];
	} else {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: request for value outside tensor.";
	}
	return values[0];
}



