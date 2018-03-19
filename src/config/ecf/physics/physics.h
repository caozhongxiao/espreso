
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICS_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICS_H_

#include "../material/material.h"

namespace espreso {

enum class PHYSICS {
	HEAT_TRANSFER_2D,
	HEAT_TRANSFER_3D,
	STRUCTURAL_MECHANICS_2D,
	STRUCTURAL_MECHANICS_3D,
	SHALLOW_WATER_2D
};

enum class DISCRETIZATION {
		FEM,
		BEM
};

struct PhysicsConfiguration: public ECFObject {

	enum class INTERPOLATION {
		LINEAR,
		QUADRATIC
	};

	size_t load_steps;

	INTERPOLATION interpolation;
	std::map<std::string, DISCRETIZATION> discretization;
	DIMENSION dimension;
	MaterialConfiguration::PHYSICAL_MODEL physical_model;

	std::map<std::string, MaterialConfiguration> materials;
	std::map<std::string, std::string> material_set;

	std::map<std::string, ECFExpression> initial_temperature, thickness;

	bool contact_interfaces;

	PhysicsConfiguration(DIMENSION dimension, MaterialConfiguration::PHYSICAL_MODEL physicalModel);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICS_H_ */
