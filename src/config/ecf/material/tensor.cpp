
#include "tensor.h"

#include "../../../basis/logging/logging.h"
#include <iostream>

using namespace espreso;

TensorConfiguration::TensorConfiguration(size_t size)
: size(size), values(size * size)
{

}

size_t TensorConfiguration::_get(size_t row, size_t column) const
{
	if (row < size && column < size) {
		return row * size + column;
	} else {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: request for value outside tensor.";
	}
	return 0;
}

const ECFExpression& TensorConfiguration::get(size_t row, size_t column) const
{
	return values[_get(row, column)];
}

ECFExpression& TensorConfiguration::get(size_t row, size_t column)
{
	return values[_get(row, column)];
}



