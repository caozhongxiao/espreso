
#ifndef SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_
#define SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_

#include "solver.h"
#include "discretization.h"
#include "../material/material.h"

namespace espreso {

struct StructuralMechanicsConfiguration: public ECFObject {

	StructuralMechanicsPhysicsSolverConfiguration physics_solver;

	std::map<std::string, std::string> material_set, initial_temperature;
	std::map<size_t, std::map<std::string, std::string> > temperature, displacement, acceleration, normal_pressure, obstacle, normal_direction;

	bool post_process;

	StructuralMechanicsConfiguration();
};

struct StructuralMechanics2DConfiguration: public StructuralMechanicsConfiguration {

	enum class ELEMENT_BEHAVIOUR {
		PLANE_STRAIN = 0,
		AXISYMMETRIC = 1,
		PLANE_STRESS = 2,
		PLANE_STRESS_WITH_THICKNESS = 3
	};

	ELEMENT_BEHAVIOUR element_behaviour;

	std::map<size_t, std::map<std::string, std::string> > thickness, angular_velocity;
	std::map<std::string, MaterialConfiguration> materials;

	StructuralMechanics2DConfiguration();
};

struct StructuralMechanics3DConfiguration: public StructuralMechanicsConfiguration {

	DISCRETIZATION discretization;
	std::map<std::string, MaterialConfiguration> materials;

	StructuralMechanics3DConfiguration();
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_ */
