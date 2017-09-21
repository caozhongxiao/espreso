
#ifndef SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_
#define SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_

#include "coordinatesystem.h"
#include "linearelasticproperties.h"
#include "thermalconductivity.h"

namespace espreso {

struct MaterialConfiguration: public ECFObject {

	enum PHYSICAL_MODEL {
		THERMAL        = 1 << 0,
		LINEAR_ELASTIC = 1 << 1
	};

	std::string name;
	std::string description;

	CoordinateSystemConfiguration coordinate_system;

	PHYSICAL_MODEL physical_model;

	ECFExpression density;
	ECFExpression heat_capacity;
	LinearElasticPropertiesConfiguration linear_elastic_properties;
	ThermalConductivityConfiguration thermal_conductivity;

	MaterialConfiguration();
	MaterialConfiguration(DIMENSION dimension, PHYSICAL_MODEL allowedPhysicalModels);

protected:
	PHYSICAL_MODEL _allowed_physical_models;
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_ */
