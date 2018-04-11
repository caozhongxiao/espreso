/*
 * hypre.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: beh01
 */

#include "hypre.h"


#include "include/HYPRE_krylov.h"
#include "include/HYPRE.h"
#include "include/HYPRE_parcsr_ls.h"


using namespace espreso;

void HYPRE::Solve(int rank, int processes, MPI_Comm communicator ) {


	 int N = 8;
	 int local_size = N/processes;

	 int ilower = local_size*rank;
	 int iupper = local_size*(rank+1) -1;

	 printf("%d: Matrix size %d, start %d end %d\n", rank, N, ilower, iupper);

	/* Create the matrix.
	      Note that this is a square matrix, so we indicate the row partition
	      size twice (since number of rows = number of cols) */
	 HYPRE_IJMatrix A;

	 HYPRE_IJMatrixCreate(communicator, ilower, iupper, ilower, iupper, &A);

	   /* Choose a parallel csr format storage (see the User's Manual) */
	 HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

	   /* Initialize before setting coefficients */
	   HYPRE_IJMatrixInitialize(A);

	   double values[3];
	   values[0] = 1;
	   values[1] = 2;
	   values[2] = 1;

	   for (int i = ilower; i <= iupper; i++)
	   {

	   /*
	HYPRE IJMatrixSetValues (HYPRE IJMatrix matrix, int nrows, int* ncols,
	const int* rows, const int* cols, const HYPRE Complex* values)*/
	    	if (i == 0) {
		    int count = 2;
		    int cols[2];
		    cols[0] = i;
	            cols[1] = i+1;
		    HYPRE_IJMatrixSetValues(A, 1, &count, &i, cols, &values[1]);
	   	} else if (i == processes*local_size -1) {
		    int count = 2;
		    int cols[2];
		    cols[0] = i-1;
	            cols[1] = i;
		    HYPRE_IJMatrixSetValues(A, 1, &count, &i, cols, values);

		} else {
		    int count = 3;
		    int cols[3];
		    cols[0] = i-1;
	            cols[1] = i;
	            cols[2] = i+1;
		    HYPRE_IJMatrixSetValues(A, 1, &count, &i, cols, values);
		}
	   }
	   HYPRE_IJMatrixAssemble(A);

	 /* Create the rhs and solution */
	   HYPRE_IJVector b;
	   HYPRE_IJVectorCreate(communicator, ilower, iupper,&b);
	   HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
	   HYPRE_IJVectorInitialize(b);

	   HYPRE_IJVector x;
	   HYPRE_IJVectorCreate(communicator, ilower, iupper,&x);
	   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
	   HYPRE_IJVectorInitialize(x);

	   double rhs_values[local_size], x_values[local_size];
	   int    rows[local_size];
	   for (int i=0; i<local_size; i++)
	   {
	      rhs_values[i] = 1.0;
	      x_values[i] = 0.0;
	      rows[i] = ilower + i;
	   }

	   HYPRE_IJVectorSetValues(b, local_size, rows, rhs_values);
	   HYPRE_IJVectorSetValues(x, local_size, rows, x_values);

	   HYPRE_IJVectorAssemble(b);

	   HYPRE_IJVectorAssemble(x);

	   HYPRE_IJMatrixPrint(A, "test.A");
	   HYPRE_IJVectorPrint(b, "test.b");

	   HYPRE_ParCSRMatrix parcsr_A;
	   HYPRE_ParVector par_b;
	   HYPRE_ParVector par_x;
	   HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
	   HYPRE_IJVectorGetObject(b, (void **) &par_b);
	   HYPRE_IJVectorGetObject(x, (void **) &par_x);

	   int num_iterations;
	   double final_res_norm;
	   HYPRE_Solver solver, precond;

	      /* Create solver */
	      HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

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
	      HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
	      HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

	      /* Run info - needed logging turned on */
	      HYPRE_PCGGetNumIterations(solver, &num_iterations);
	      HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
	      if (rank == 0)
	      {
	         printf("\n");
	         printf("Iterations = %d\n", num_iterations);
	         printf("Final Relative Residual Norm = %e\n", final_res_norm);
	         printf("\n");
	      }

	      /* Destroy solver and preconditioner */
	      HYPRE_ParCSRPCGDestroy(solver);
	      HYPRE_BoomerAMGDestroy(precond);

	   HYPRE_IJVectorPrint(x, "test.x");

	   HYPRE_IJMatrixDestroy(A);
	   HYPRE_IJVectorDestroy(b);
	   HYPRE_IJVectorDestroy(x);
}
