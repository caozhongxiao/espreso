
#include "hypre.h"

#include "../../basis/logging/logging.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"

#include "include/HYPRE_krylov.h"
#include "include/HYPRE.h"
#include "include/HYPRE_parcsr_ls.h"

#include <vector>
#include <numeric>
#include "../../config/ecf/linearsolver/hypre/hypre.h"

using namespace espreso;

HypreData::HypreData(MPI_Comm &comm, eslocal nrows)
: _comm(comm), _roffset(nrows), _nrows(nrows), _finalized(false)
{
	Communication::exscan(_roffset);
	HYPRE_IJMatrixCreate(comm, _roffset + 1, _roffset + _nrows, _roffset + 1, _roffset + _nrows, &_K);
	HYPRE_IJVectorCreate(comm, _roffset + 1, _roffset + _nrows, &_f);
	HYPRE_IJVectorCreate(comm, _roffset + 1, _roffset + _nrows, &_x);

	HYPRE_IJMatrixSetObjectType(_K, HYPRE_PARCSR);
	HYPRE_IJVectorSetObjectType(_f, HYPRE_PARCSR);
	HYPRE_IJVectorSetObjectType(_x, HYPRE_PARCSR);

	HYPRE_IJMatrixInitialize(_K);
	HYPRE_IJVectorInitialize(_f);
	HYPRE_IJVectorInitialize(_x);
}

void HypreData::insertCSR(eslocal nrows, eslocal offset, eslocal *rowPrts, eslocal *colIndices, double *values, double *rhsValues)
{
	std::vector<eslocal> ncols, rows;
	ncols.reserve(nrows);
	rows.reserve(nrows);
	for (eslocal r = 0; r < nrows; r++) {
		ncols.push_back(rowPrts[r + 1] - rowPrts[r]);
		rows.push_back(_roffset + offset + r + 1);
	}
	std::vector<double> x(nrows);

	HYPRE_IJMatrixSetValues(_K, nrows, ncols.data(), rows.data(), colIndices, values);
	HYPRE_IJVectorSetValues(_f, nrows, rows.data(), rhsValues);
	HYPRE_IJVectorSetValues(_x, nrows, rows.data(), x.data());

	_finalized = false;
}

void HypreData::insertIJV(eslocal nrows, eslocal offset, eslocal size, eslocal *rowIndices, eslocal *colIndices, double *values, double *rhsValues)
{
	std::vector<eslocal> ncols, rows;
	ncols.reserve(nrows);
	rows.reserve(nrows);
	for (eslocal i = 0; i < size; i++) {
		if (ncols.size() == 0 || ncols.back() != colIndices[i]) {
			ncols.push_back(1);
		} else {
			++ncols.back();
		}
		rows.push_back(_roffset + offset + rowIndices[i]);
	}
	std::vector<double> x(nrows);

	HYPRE_IJMatrixSetValues(_K, nrows, ncols.data(), rows.data(), colIndices, values);
	HYPRE_IJVectorSetValues(_f, nrows, rows.data(), rhsValues);
	HYPRE_IJVectorSetValues(_x, nrows, rows.data(), x.data());
	_finalized = false;
}

void HypreData::finalizePattern()
{
	HYPRE_IJMatrixAssemble(_K);
	HYPRE_IJVectorAssemble(_f);
	HYPRE_IJVectorAssemble(_x);
	_finalized = true;
}

HypreData::~HypreData()
{
	HYPRE_IJMatrixDestroy(_K);
	HYPRE_IJVectorDestroy(_f);
	HYPRE_IJVectorDestroy(_x);
}

static void setBoomerAMG(HYPRE_Solver &boomerAMG, const HYPREBoomerAMGConfiguration &configuration)
{
	HYPRE_BoomerAMGCreate(&boomerAMG);

	HYPRE_BoomerAMGSetPrintLevel(boomerAMG, 0);
	HYPRE_BoomerAMGSetCoarsenType(boomerAMG, 6);
	HYPRE_BoomerAMGSetOldDefault(boomerAMG);
	HYPRE_BoomerAMGSetRelaxType(boomerAMG, 6);
	HYPRE_BoomerAMGSetNumSweeps(boomerAMG, 1);
	HYPRE_BoomerAMGSetTol(boomerAMG, 0.0);
	HYPRE_BoomerAMGSetMaxIter(boomerAMG, 1);
}

void HYPRE::solve(const HypreConfiguration &configuration, HypreData &data, eslocal nrows, double *solution)
{
	if (!data._finalized) {
		data.finalizePattern();
	}

	HYPRE_ParCSRMatrix K;
	HYPRE_ParVector f, x;
	HYPRE_IJMatrixGetObject(data._K, (void**) &K);
	HYPRE_IJVectorGetObject(data._f, (void**) &f);
	HYPRE_IJVectorGetObject(data._x, (void**) &x);

	HYPRE_Solver solver;
	HYPRE_Solver preconditioner;
	switch (configuration.solver_type) {
	case HypreConfiguration::SOLVER_TYPE::BoomerAMG:
		setBoomerAMG(solver, configuration.boomeramg);
		break;
	case HypreConfiguration::SOLVER_TYPE::PCG:
		HYPRE_ParCSRPCGCreate(data._comm, &solver);

//		HYPRE_PCGSetMaxIter(solver, configuration.max_iterations);
		HYPRE_PCGSetTol(solver, configuration.pcg.relative_conv_tol);
		HYPRE_PCGSetTwoNorm(solver, 1);
		HYPRE_PCGSetPrintLevel(solver, 0);
		HYPRE_PCGSetLogging(solver, 0);

		switch (configuration.pcg.preconditioner) {
		case HYPREPCGConfiguration::PRECONDITIONER::BoomerAMG:
			setBoomerAMG(preconditioner, configuration.pcg.boomeramg);

			HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
			break;
		case HYPREPCGConfiguration::PRECONDITIONER::Parasalis:
		default:
			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: not implemented interface to the required solver.";
		}

		break;
	default:
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: not implemented interface to the required solver.";
	}

//	HYPRE_IJMatrixPrint(data._K, "test.K");
//	HYPRE_IJVectorPrint(data._f, "test.f");

	HYPRE_ParCSRPCGSetup(solver, K, f, x);
	HYPRE_ParCSRPCGSolve(solver, K, f, x);

//	HYPRE_IJVectorPrint(data._x, "test.x");

	eslocal iterations;
	double norm;
	HYPRE_PCGGetNumIterations(solver, &iterations);
	HYPRE_PCGGetFinalRelativeResidualNorm(solver, &norm);

	ESINFO(CONVERGENCE) << "Final Relative Residual Norm " << norm << " in " << iterations << " iteration.";

	std::vector<eslocal> rows(nrows);
	std::iota(rows.begin(), rows.end(), data._roffset + 1);
	HYPRE_IJVectorGetValues(data._x, data._nrows, rows.data(), solution);
}
