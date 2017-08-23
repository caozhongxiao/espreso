
#ifndef SRC_CONFIG_ECF_PHYSICS_SOLVER_H_
#define SRC_CONFIG_ECF_PHYSICS_SOLVER_H_

#include "loadsteps.h"

namespace espreso {

struct PhysicsSolverConfiguration: public ECFObject {

	enum class INTERPOLATION {
		LINEAR,
		QUADRATIC
	};

	INTERPOLATION interpolation;
	size_t load_steps;

	PhysicsSolverConfiguration();
};

struct AdvectionDiffusionPhysicsSolverConfiguration: PhysicsSolverConfiguration {
	std::map<size_t, AdvectionDiffusionLoadStepsConfiguration> load_steps_settings;

	AdvectionDiffusionPhysicsSolverConfiguration();
};

struct StructuralMechanicsPhysicsSolverConfiguration: PhysicsSolverConfiguration {
	std::map<size_t, StructuralMechanicsLoadStepsConfiguration> load_steps_settings;

	StructuralMechanicsPhysicsSolverConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_PHYSICS_SOLVER_H_ */
