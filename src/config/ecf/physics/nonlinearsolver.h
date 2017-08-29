
#ifndef SRC_CONFIG_ECF_PHYSICS_NONLINEARSOLVER_H_
#define SRC_CONFIG_ECF_PHYSICS_NONLINEARSOLVER_H_

#include "../../configuration.h"

namespace espreso {

struct NonLinearSolverConfiguration: public ECFObject {

	enum class METHOD {
		NEWTON_RHAPSON,
		MODIFIED_NEWTON_RHAPSON
	};

	enum class STEPPINGG {
		TRUE,
		FALSE,
		AUTO
	};

	METHOD method;
	STEPPINGG stepping;

	size_t max_iterations, substeps;
	bool line_search, tangent_matrix_correction, adaptive_precision;

	bool check_first_residual, check_second_residual;
	double requested_first_residual, requested_second_residual;

	double r_tol, c_fact;

	std::string solution_name, residual_name;
	NonLinearSolverConfiguration(const std::string &solution_name, const std::string &residual_name);
};

struct AdvectionDiffusionNonLinearSolverConfiguration: public NonLinearSolverConfiguration {

	AdvectionDiffusionNonLinearSolverConfiguration(): NonLinearSolverConfiguration("temperature", "heat") {}
};

struct StructuralMechanicsNonLinearSolverConfiguration: public NonLinearSolverConfiguration {

	StructuralMechanicsNonLinearSolverConfiguration(): NonLinearSolverConfiguration("displacement", "forces") {}
};

}



#endif /* SRC_CONFIG_ECF_PHYSICS_NONLINEARSOLVER_H_ */
