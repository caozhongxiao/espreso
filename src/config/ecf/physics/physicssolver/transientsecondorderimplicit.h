
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTSECONDORDERIMPLICIT_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTSECONDORDERIMPLICIT_H_

#include <config/ecf/physics/physicssolver/transientfirstorderimplicit.h>

namespace espreso {

struct TransientSecondOrderImplicitConfiguration: public ECFDescription {

	enum class METHOD {
		NEWMARK
	};

	enum class MASS_MATRIX_TYPE {
		CONSISTENT,
		DIAGONAL,
		HRZDIAGONAL
	};

	enum class DUMPING {
		NONE,
		DIRECT,
		DUMPING_RATIO
	};

	METHOD method;
	AutoTimeSteppingConfiguration auto_time_stepping;
	double alpha, delta, numerical_dumping, time_step;

	DUMPING dumping;
	double mass_dumping, stiffness_dumping;
	double dumping_ratio, frequency;

	MASS_MATRIX_TYPE mass_matrix_type;

	TransientSecondOrderImplicitConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_TRANSIENTSECONDORDERIMPLICIT_H_ */
