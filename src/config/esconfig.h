
#ifndef ESCONFIG_H_
#define ESCONFIG_H_

#include <cstdlib>
#include <string>
#include <vector>
#include <map>

namespace espreso {

class Parameter;

namespace config {

extern std::vector<Parameter> parameters;


namespace env {
	extern int MPIrank;
	extern int MPIsize;

	extern size_t MKL_NUM_THREADS;
	extern size_t OMP_NUM_THREADS;
	extern size_t SOLVER_NUM_THREADS;
	extern size_t PAR_NUM_THREADS;
	extern size_t CILK_NWORKERS;

	extern std::string executable;
	extern std::string configurationFile;
};

namespace mesh {
	enum class INPUTalternative {
		/// Ansys generated by Matsol library
		MATSOL = 0,
		/// Ansys Workbench format
		WORKBENCH = 1,
		/// OpenFOAM format
		OPENFOAM = 2,
		/// ESPRESO binary format
		ESDATA = 3,
		/// ESPRESO internal problem generator
		GENERATOR = 4
	};
	/// The type of an input problem.
	/**
	 * It specifies the format of an input.
	 * Formats that support only few MPI processes (e.g. Ansys) can be converted
	 * to ESPRESO format by 'decomposer'
	 */
	extern INPUTalternative INPUT;


	/// A path to an input problem
	/**
	 * A path to an input file in case of Workbench or ESPRESO Problem Generator.
	 * In case of other format, a path to the root director of an input probem.
	 */
	extern std::string PATH;

	/// The number of sub-domains in each cluster.
	/**
	 * The number of sub-domains should be set appropriate to a problem size.
	 * Big sub-domains increase time to solve them. Small sub-domains increase
	 * the size of the coarse problem.
	 *
	 * The shared memory parallelization is mainly through sub-domains. Hence,
	 * the minimal number should be kept higher that the number of threads.
	 */
	extern size_t SUBDOMAINS;

	/// The number of fix points in each sub-domain.
	/**
	 * Fix points are nodes used in the regularization process of stiffness matrix.
	 * It is not recommended to change the default value.
	 */
	extern size_t FIX_POINTS;

	/// The number of selected points on edges and faces.
	/**
	 * TODO: Alex - co presne delaji cornery
	 */
	extern size_t CORNERS;

	/// All vertex points are marked as corner points.
	/**
	 * The 'vertex' is a end point of the 'edge'. If the 'edge' has not any end
	 * ('edge' is a circle), 4 uniformly distributed points are marked.
	 */
	extern bool VERTEX_CORNERS;

	/// Uniformly distributed points on all edges are marked as corner points.
	/**
	 * The 'edges' are border lines of the 'faces'.
	 * The number of corner points on each edge is determined by parameter 'CORNERS'.
	 */
	extern bool EDGE_CORNERS;

	/// Uniformly distributed points on all faces are marked as corner points.
	/**
	 * The 'face' is a common part between two sub-domains.
	 * The number of corner points on each face is determined by parameter 'CORNERS'.
	 */
	extern bool FACE_CORNERS;

	/// Values at the stiffness matrix corresponding to points on ‘edges’ are averaged.
	/**
	 * ESPRESO contains only a premature state of the averaging. Always keep the default value!
	 */
	extern bool AVERAGE_EDGES;

	/// Values at the stiffness matrix corresponding to points on ‘faces’ are averaged.
	/**
	 * ESPRESO contains only a premature state of the averaging. Always keep the default value!
	 */
	extern bool AVERAGE_FACES;
};

namespace output {

	enum class OUTPUT_FORMATAlternatives {
		VTK_LEGACY_FORMAT = 0,
		VTK_BINARY_FORMAT = 1,
		VTK_MULTIBLOCK_FORMAT = 2,
		ENSIGHT_FORMAT = 3
	};
	/// Format of output data
	extern OUTPUT_FORMATAlternatives OUTPUT_FORMAT;

	/// All results are compressed by 'z' library
	extern bool OUTPUT_COMPRESSION;

	/// Mesh is decimated by this ratio
	extern double OUTPUT_DECIMATION;

	/// Save nodes with Dirichlet condition to VTK files.
	extern bool SAVE_PROPERTIES;

	/// Save gluing of sub-domains and clusters
	extern bool SAVE_GLUING;

	/// Save result computed by the solver to VTK files.
	extern bool SAVE_RESULTS;


	/// All sub-domains are shrunk by this ratio.
	extern double SUBDOMAINS_SHRINK_RATIO;

	/// All clusters are shrunk by this ratio.
	extern double CLUSTERS_SHRINK_RATIO;
};

namespace assembler {
	enum class DISCRETIZATIONalternative {
		/// Finite Element Method
		FEM = 0,
		/// Boundary Element Method
		BEM = 1
	};
	/// Discretization of an example.
	/**
	 * Stiffness matrices are computed based on the discretization:
	 * - Finite Element Method is used
	 * - Boundary Element Method is used
	 */
	extern DISCRETIZATIONalternative DISCRETIZATION;

	enum class DOFS_ORDERalternative {
		/// Group elements - x1, y1, z1, x2, y2, z2, ....
		GROUP_ELEMENTS = 0,

		/// Group elements - x1, x2, ..., y1, y2, ..., z1, z2, ....
		GROUP_DOFS = 1
	};

	extern DOFS_ORDERalternative DOFS_ORDER;
};

namespace solver {
	/// ESPRESO checks the norm of the solution if NORM is not zero
	extern double NORM;

	/// The solver requested precision.
	extern double EPSILON;

	/// Maximum iterations for the solver.
	extern size_t ITERATIONS;


	enum class FETI_METHODalternative {
		/// Total FETI
		TOTAL_FETI = 0,
		/// Hybrid Total FETI
		HYBRID_FETI = 1,
		/// Multi-grid Hypre interface
		HYPRE = 2
	};
	/// A variant of FETI method used by the solver
	extern FETI_METHODalternative FETI_METHOD;


	enum class PRECONDITIONERalternative {
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
	/// Used preconditioner
	extern PRECONDITIONERalternative PRECONDITIONER;


	enum class REGULARIZATIONalternative {
		/// Fix points
		FIX_POINTS = 0,
		/// Randomly found null pivots of stiffness matrix
		NULL_PIVOTS = 1
	};

	/// A type of regularization of stiffness matrix.
	/**
	 * In the case of a singular stiffness matrix, regularization has to be applied.
	 * When an example is loaded with mesh, it is possible to use regularization
	 * from fix points. It is faster then regularization from null pivots.
	 * However, null pivots are more general and usable even without mesh
	 * (e.g. when API is used).
	 */
	extern REGULARIZATIONalternative REGULARIZATION;

	/// Use redundant Lagrange multipliers.
	/**
	 * In the case of Hybrid Total FETI, ESPRESO compose gluing matrix for
	 * each cluster. If this option is on, the multipliers from the cluster
	 * gluing matrix will be also in global gluing matrix.
	 */
	extern bool REDUNDANT_LAGRANGE;

	enum class B0_TYPEalternative {
		/// Gluing based on corners
		CORNERS = 0,
		/// Gluing based on kernels of faces
		KERNELS = 1,
		/// Both corners and kernels
		COMBINED = 2
	};
	/// Type of cluster gluing matrix.
	/**
	 * Sub-domains in each cluster have to be glued together. Gluing can be based
	 * on corners (random nodes at interfaces between sub-domains) or kernels of
	 * faces between sub-domains.
	 */
	extern B0_TYPEalternative B0_TYPE;


	/// Schur complement will be used.
	extern bool USE_SCHUR_COMPLEMENT;

	enum class SCHUR_COMPLEMENT_PRECalternative {
		/// Double precision
		DOUBLE = 0,
		/// Single precision
		SINGLE = 1
	};
	/// Precision of Schur complement
	extern SCHUR_COMPLEMENT_PRECalternative SCHUR_COMPLEMENT_PREC;

	enum class SCHUR_COMPLEMENT_TYPEalternative {
		/// A full matrix is stored
		GENERAL = 0,
		/// Store only triangle
		SYMMETRIC = 1
	};
	extern SCHUR_COMPLEMENT_TYPEalternative SCHUR_COMPLEMENT_TYPE;

	/// Combine usage of SC for Accelerator and Sparse Direct Solver for CPU
	extern bool COMBINE_SC_AND_SPDS;

	/// Keep factors between iterations
	extern bool KEEP_FACTORS;


	enum class CGSOLVERalternative {
		STANDARD = 0,
		PIPELINED = 1,
		FULL_ORTOGONAL = 2,
		GMRES = 3,
		BICGSTAB = 4,
		QPCE = 5
	};

	/// A type of conjugate gradient solver
	extern CGSOLVERalternative CGSOLVER;


	enum class KSOLVERalternative {
		/// A direct solver with double precision
		DIRECT_DP = 0,
		/// An iterative solver
		ITERATIVE = 1,
		/// A direct solver with single precision
		DIRECT_SP = 2,
		/// A direct solver with mixed precision
		DIRECT_MP = 3
	};
	/// A type of stiffness matrix solver
	extern KSOLVERalternative KSOLVER;

	/// Number of reiteration steps for single precision direct solver
	extern size_t   KSOLVER_SP_STEPS;
	/// Iterative norm single precision direct solver
	extern double   KSOLVER_SP_NORM;


	enum class F0SOLVERalternative {
		/// The same precision as K solver
		K_PRECISION = 0,
		/// Double precision
		DOUBLE = 1
	};
	/// A type of F0 solver
	extern F0SOLVERalternative F0SOLVER;

	enum class SASOLVERalternative {
		CPU_DENSE = 0,
		ACC_DENSE = 1,
		CPU_SPARSE = 2
	};
	/// A type of S alfa solver
	extern SASOLVERalternative SASOLVER;

	/// Number of used MIC accelerators
	extern size_t N_MICS;

	/// The number of time steps for transient problems
	extern size_t TIME_STEPS;
};

namespace hypre {
	enum class SOLVERalternative { 
		CG=0,
		GMRES = 1,
		FGMRES = 2,
		BICGS = 3,
		BICGSTAB = 4,
		TFQMR = 5,
		SYMQMR = 6,
		SUPERLU = 7,
		SUPERLUX = 8
	}; 
	
	enum class PRECONDITIONERalternative {
		DIAGONAL = 0,
		PILUT = 1,
		EUCLID = 2,
		PARASAILS = 3,
		BOOMERAMG = 4,
		POLY = 5,
		MLI = 6
	};
	
	extern SOLVERalternative HYPRE_SOLVER;
	extern PRECONDITIONERalternative HYPRE_PRECONDITIONER;
};

namespace info {
	/// An output directory for files generated by debug run.
	/**
	 * When printMatrices is set, ESPRESO saves matrices to this directory.
	 * Matrices are divided into directories according to MPI processes.
	 * All matrices are printed in format {matrix}{sub-domain}.txt.
	 */
	extern std::string OUTPUT;

	/// ESPRESO verbose level
	/**
	 * ESPRESO print an info based on verbose level. The default value is 0
	 * and maximal value is 3. Higher value does not add any more information.
	 */
	extern size_t VERBOSE_LEVEL;

	/// ESPRESO testing level
	/**
	 * ESPRESO contains an infrastructure for on-the-fly testing.
	 * Trivially testable values are always tested. Increasing the testing
	 * level turn time dependent testing on.
	 */
	extern size_t TESTING_LEVEL;

	/// ESPRESO measure level
	/**
	 * Measure the time for computing particular functions. In default, all
	 * measurements are disabled.
	 */
	extern size_t MEASURE_LEVEL;

	/// The main ESPRESO debug feature. All matrices are printed.
	/**
	 * ESPRESO print all matrices to 'output' directory. The FETI computation
	 * can be reproduces from these matrices by run 'espreso.py' in python
	 * directory.
	 */
	extern bool PRINT_MATRICES;
};

}

}


#endif /* ESCONFIG_H_ */
