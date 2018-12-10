
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_FETI_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_FETI_H_

#include "../../configuration.h"

namespace espreso {

enum class FETI_METHOD {
	/// Total FETI
	TOTAL_FETI = 0,
	/// Hybrid Total FETI
	HYBRID_FETI = 1,
};

enum class FETI_ITERATIVE_SOLVER {
	/// Projected conjugate gradients
	PCG = 0,
	/// Pipelined PCG
	pipePCG = 1,
	/// Full orthogonal PCG
	orthogonalPCG = 2,
	/// GMRES
	GMRES = 3,
	/// BICGSTAB
	BICGSTAB = 4,
	/// QPCE
	QPCE = 5,
	/// FETI Geneo with full ortogonalization CG
	orthogonalPCG_CP = 6,
	/// FETI Geneo with regular CG
	PCG_CP = 7
};

enum class FETI_PRECONDITIONER {
	/// No preconditioner is used
	NONE = 0,
	/// Lumped preconditioner     S = K_ss
	LUMPED = 1,
	/// Weight function
	WEIGHT_FUNCTION = 2,
	/// Dirichlet preconditioner  S = K_ss - K_sr * inv(K_rr) * K_sr
	DIRICHLET = 3,
	/// simplified Dirichlet      S = K_ss - K_sr * 1/diag(K_rr) * K_sr
	SUPER_DIRICHLET = 4,
	/// Lubos's preconditioner
	MAGIC = 5
};

enum class FETI_CONJ_PROJECTOR {
	/// No conj projector
	NONE = 0,
	/// Randomly found null pivots of stiffness matrix
	GENEO,
	/// Conjugate projector for transient problems
	CONJ_R,
	CONJ_K
};

enum class FETI_REGULARIZATION {
	/// Based on a physics
	ANALYTIC = 0,
	/// Randomly found null pivots of stiffness matrix
	ALGEBRAIC = 1
};

enum class FETI_B0_TYPE {
	/// Gluing based on corners
	CORNERS = 0,
	/// Gluing based on kernels of faces
	KERNELS = 1,
};

enum class FETI_FLOAT_PRECISION {
	/// Double precision
	DOUBLE = 0,
	/// Single precision
	SINGLE = 1
};

enum class FETI_KSOLVER {
	/// A direct solver with double precision
	DIRECT_DP = 0,
	/// An iterative solver
	ITERATIVE = 1,
	/// A direct solver with single precision
	DIRECT_SP = 2,
	/// A direct solver with mixed precision
	DIRECT_MP = 3
};

enum class FETI_F0SOLVER_PRECISION {
	/// The same precision as K solver
	K_PRECISION = 0,
	/// Double precision
	DOUBLE = 1
};

enum class FETI_SASOLVER {
	/// Dense on CPU
	CPU_DENSE = 0,
	/// Dense on accelerator
	ACC_DENSE = 1,
	/// Sparse on CPU
	CPU_SPARSE = 2
};

enum class FETI_MATRIX_STORAGE {
	/// A full matrix is stored
	GENERAL = 0,
	/// Store only triangle
	SYMMETRIC = 1
};

struct FETISolverConfiguration: public ECFObject {

	double precision;
	size_t max_iterations;

	FETI_METHOD method;
	FETI_ITERATIVE_SOLVER iterative_solver;
	FETI_PRECONDITIONER preconditioner;
	FETI_REGULARIZATION regularization;
	FETI_CONJ_PROJECTOR conjugate_projector;

	size_t geneo_size, restart_iteration, num_restart;

	bool orthogonal_K_kernels;
	bool redundant_lagrange, scaling;

	FETI_B0_TYPE B0_type;

	bool use_schur_complement;

	FETI_FLOAT_PRECISION schur_precision;
	FETI_KSOLVER Ksolver;
	size_t Ksolver_max_iterations;
	double Ksolver_precision;
	FETI_F0SOLVER_PRECISION F0_precision;
	FETI_SASOLVER SAsolver;
	FETI_MATRIX_STORAGE schur_type;

	bool mp_pseudoinverse, combine_sc_and_spds, keep_factors;

	size_t sc_size, n_mics;
	bool load_balancing, load_balancing_preconditioner;

	FETISolverConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_FETI_H_ */
