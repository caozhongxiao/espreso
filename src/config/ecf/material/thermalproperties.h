
#ifndef SRC_CONFIG_ECF_MATERIAL_THERMALPROPERTIES_H_
#define SRC_CONFIG_ECF_MATERIAL_THERMALPROPERTIES_H_

#include "tensor.h"

namespace espreso {

struct ThermalPropertiesConfiguration: public ECFObject {

	enum class MODEL {
		ISOTROPIC,
		DIAGONAL,
		SYMMETRIC,
		ANISOTROPIC,
	};

	MODEL model;
	bool is3D;

	TensorConfiguration thermal_conductivity;

	ThermalPropertiesConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_THERMALPROPERTIES_H_ */
