
#ifndef SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_
#define SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_

#include "coordinatesystem.h"
#include "linearelasticproperties.h"
#include "thermalproperties.h"

namespace espreso {

struct MaterialConfiguration: public ECFObject {

	enum PHYSICAL_MODEL {
		THERMAL        = 1 << 0,
		LINEAR_ELASTIC = 1 << 1
	};

	PHYSICAL_MODEL physical_model;

	ECFExpression density;
	ECFExpression heat_capacity;
	LinearElasticPropertiesConfiguration linear_elastic_properties;
	ThermalPropertiesConfiguration thermal_properties;

	MaterialConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_ */
