
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTSECONDORDERIMPLICIT_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTSECONDORDERIMPLICIT_H_

#include "transientfirstorderimplicit.h"
#include "config/holders/expression.h"

namespace espreso {

struct TransientSecondOrderImplicitDirectDampingConfiguration: public ECFDescription {

	ECFExpression mass, stiffness;

	TransientSecondOrderImplicitDirectDampingConfiguration();
};

struct TransientSecondOrderImplicitRatioDampingConfiguration: public ECFDescription {

	ECFExpression ratio, frequency;

	TransientSecondOrderImplicitRatioDampingConfiguration();
};

struct TransientSecondOrderImplicitSolverConfiguration: public ECFDescription {

	enum class METHOD {
		NEWMARK
	};

	enum class MASS_MATRIX_TYPE {
		CONSISTENT,
		DIAGONAL,
		HRZDIAGONAL
	};

	enum class DAMPING {
		NONE,
		DIRECT,
		DAMPING_RATIO
	};

	METHOD method;
	AutoTimeSteppingConfiguration auto_time_stepping;
	double alpha, delta, numerical_damping, time_step;

	DAMPING damping;

	TransientSecondOrderImplicitDirectDampingConfiguration direct_damping;
	TransientSecondOrderImplicitRatioDampingConfiguration ratio_damping;

	MASS_MATRIX_TYPE mass_matrix_type;

	TransientSecondOrderImplicitSolverConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTSECONDORDERIMPLICIT_H_ */
