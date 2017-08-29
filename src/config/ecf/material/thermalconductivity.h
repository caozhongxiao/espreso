
#ifndef SRC_CONFIG_ECF_MATERIAL_THERMALCONDUCTIVITY_H_
#define SRC_CONFIG_ECF_MATERIAL_THERMALCONDUCTIVITY_H_

#include "tensor.h"

namespace espreso {

struct ThermalConductivityConfiguration: public ECFObject {

	enum class MODEL {
		ISOTROPIC,
		DIAGONAL,
		SYMMETRIC,
		ANISOTROPIC,
	};

	MODEL model;
	bool is3D;

	TensorConfiguration values;

	ThermalConductivityConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_THERMALCONDUCTIVITY_H_ */
