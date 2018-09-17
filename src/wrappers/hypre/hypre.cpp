
#include "hypre.h"

#include "include/HYPRE_krylov.h"
#include "include/HYPRE.h"
#include "include/HYPRE_parcsr_ls.h"

using namespace espreso;

void HYPRE::Solve(MultigridConfiguration &configuration, int rank,
		int processes, MPI_Comm communicator,
		esglobal start_row, eslocal num_rows, std::vector<HypreRegion> values,
		double* f_data, double* result_data) {

	esglobal end_row = start_row + num_rows-1;


	/* Create the matrix.
	 Note that this is a square matrix, so we indicate the row partition
	 size twice (since number of rows = number of cols) */
	HYPRE_IJMatrix K;

	HYPRE_IJMatrixCreate(communicator, start_row, end_row, start_row, end_row, &K);

	/* Choose a parallel csr format storage (see the User's Manual) */
	HYPRE_IJMatrixSetObjectType(K, HYPRE_PARCSR);

	/* Initialize before setting coefficients */
	HYPRE_IJMatrixInitialize(K);

	/*HYPRE IJMatrixSetValues (HYPRE IJMatrix matrix, int nrows, int* ncols,
	 const int* rows, const int* cols, const HYPRE Complex* values)*/

	for (auto region = values.begin(); region != values.end(); region++) {
		HYPRE_IJMatrixSetValues(K, region->nrows, region->ncols.data(),
				region->rows.data(), region->cols, region->values);
	}
	HYPRE_IJMatrixAssemble(K);

	/* Create the rhs and solution */
	HYPRE_IJVector f;
	HYPRE_IJVectorCreate(communicator, start_row, end_row, &f);
	HYPRE_IJVectorSetObjectType(f, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(f);

	HYPRE_IJVector result;
	HYPRE_IJVectorCreate(communicator, start_row, end_row, &result);
	HYPRE_IJVectorSetObjectType(result, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(result);

	int rows[num_rows];
	for (int i = 0; i < num_rows; i++) {
		result_data[i] = 0.0;
		rows[i] = start_row + i;
	}

	HYPRE_IJVectorSetValues(f, num_rows, rows, f_data);
	HYPRE_IJVectorSetValues(result, num_rows, rows, result_data);

	HYPRE_IJVectorAssemble(f);
	HYPRE_IJVectorAssemble(result);

//	HYPRE_IJMatrixPrint(K, "test.K");
//	HYPRE_IJVectorPrint(f, "test.f");
//	HYPRE_IJVectorPrint(result, "test.result");

	HYPRE_ParCSRMatrix parcsr_K;
	HYPRE_ParVector par_f;
	HYPRE_ParVector par_result;
	HYPRE_IJMatrixGetObject(K, (void**) &parcsr_K);
	HYPRE_IJVectorGetObject(f, (void **) &par_f);
	HYPRE_IJVectorGetObject(result, (void **) &par_result);

	int num_iterations;
	double final_res_norm;
	HYPRE_Solver solver, precond;

	/* Create solver */
	HYPRE_ParCSRPCGCreate(communicator, &solver);

	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
	HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
	HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
	HYPRE_PCGSetPrintLevel(solver, 2); /* print solve info */
	HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

	/* Now set up the AMG preconditioner and specify any parameters */
	HYPRE_BoomerAMGCreate(&precond);
	HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
	HYPRE_BoomerAMGSetCoarsenType(precond, 6);
	HYPRE_BoomerAMGSetOldDefault(precond);
	HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
	HYPRE_BoomerAMGSetNumSweeps(precond, 1);
	HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
	HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

	/* Set the PCG preconditioner */
	HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

	/* Now setup and solve! */
	HYPRE_ParCSRPCGSetup(solver, parcsr_K, par_f, par_result);
	HYPRE_ParCSRPCGSolve(solver, parcsr_K, par_f, par_result);

	/* Run info - needed logging turned on */
	HYPRE_PCGGetNumIterations(solver, &num_iterations);
	HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

	printf("\n");
	printf("Iterations = %d\n", num_iterations);
	printf("Final Relative Residual Norm = %e\n", final_res_norm);
	printf("\n");

	/* Destroy solver and preconditioner */
	HYPRE_ParCSRPCGDestroy(solver);
	HYPRE_BoomerAMGDestroy(precond);

//	HYPRE_IJVectorPrint(result, "test.result");

	HYPRE_IJVectorGetValues(result, num_rows, rows, result_data);

	HYPRE_IJMatrixDestroy(K);
	HYPRE_IJVectorDestroy(f);
	HYPRE_IJVectorDestroy(result);
}
