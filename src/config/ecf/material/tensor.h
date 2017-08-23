
#ifndef SRC_CONFIG_ECF_MATERIAL_TENSOR_H_
#define SRC_CONFIG_ECF_MATERIAL_TENSOR_H_

#include "../../configuration.h"

namespace espreso {

struct TensorConfiguration: public ECFObject {

	size_t size;
	std::vector<ECFExpression> values;

	TensorConfiguration(size_t size);

	ECFExpression& value(size_t row, size_t column);
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_TENSOR_H_ */
