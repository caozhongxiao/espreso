
#ifndef SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_
#define SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_

#include "physics.h"
#include "physicssolver/loadstep.h"

namespace espreso {

struct StructuralMechanicsGlobalSettings {

	enum class ELEMENT_BEHAVIOUR {
		PLANE_STRAIN = 0,
		AXISYMMETRIC = 1,
		PLANE_STRESS = 2,
		PLANE_STRESS_WITH_THICKNESS = 3
	};

	ELEMENT_BEHAVIOUR element_behaviour;

	std::map<std::string, ECFExpression> initial_temperature, thickness;
};

struct StructuralMechanicsLoadStepConfiguration: public StructuralMechanicsLoadStepSolverConfiguration {

	bool large_displacement;

	std::map<std::string, ECFExpression> temperature, normal_pressure;
	std::map<std::string, ECFExpressionVector> angular_velocity, acceleration, normal_direction, obstacle;
	std::map<std::string, ECFExpressionOptionalVector> displacement;

	StructuralMechanicsLoadStepConfiguration(DIMENSION *D);
};

struct StructuralMechanicsOutputSettings {

	bool displacement;

	void basic() {
		displacement = true;
	}
	void all() {
		displacement = true;
	}

	StructuralMechanicsOutputSettings() { basic(); }
};

struct StructuralMechanicsConfiguration: public PhysicsConfiguration, public StructuralMechanicsGlobalSettings {

	std::map<size_t, StructuralMechanicsLoadStepConfiguration> load_steps_settings;

	StructuralMechanicsConfiguration(DIMENSION D);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_ */
