
#ifndef SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_
#define SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_

#include "coordinatesystem.h"
#include "linearelasticproperties.h"
#include "thermalconductivity.h"

namespace espreso {

struct MaterialBaseConfiguration: public ECFObject {

	enum PHYSICAL_MODEL {
		THERMAL        = 1 << 0,
		LINEAR_ELASTIC = 1 << 1
	};

	CoordinateSystemConfiguration coordinate_system;

	ECFExpression density;
	ECFExpression heat_capacity;
	LinearElasticPropertiesConfiguration linear_elastic_properties;
	ThermalConductivityConfiguration thermal_conductivity;

	MaterialBaseConfiguration();
	MaterialBaseConfiguration(bool *phase_change, PHYSICAL_MODEL *physicalModel);
	MaterialBaseConfiguration(bool *phase_change, PHYSICAL_MODEL *physicalModel, DIMENSION *dimension);

protected:
	bool *_phase_change;
	PHYSICAL_MODEL *_physical_model;
};

struct MaterialConfiguration: public MaterialBaseConfiguration {

	std::string name;
	std::string description;

	DIMENSION dimension;
	PHYSICAL_MODEL physical_model;

	bool phase_change;
	size_t smooth_step_order;
	double latent_heat, transition_interval, phase_change_temperature;

	std::map<size_t, MaterialBaseConfiguration> phases;

	MaterialConfiguration();
	MaterialConfiguration(DIMENSION dimension, PHYSICAL_MODEL allowedPhysicalModels);

	// allow assign operator (needed because phases is std::map of ECFOBjects)
	MaterialConfiguration& operator=(const MaterialConfiguration &other);

protected:
	PHYSICAL_MODEL _allowed_physical_models;
};

}



#endif /* SRC_CONFIG_ECF_MATERIAL_MATERIAL_H_ */
