
#include "hypre.h"

#include "../../basis/logging/logging.h"
#include "../../basis/utilities/communication.h"
#include "../../basis/utilities/utils.h"

#include "../../globals/run.h"

#include "../../config/ecf/root.h"

#include "include/HYPRE_krylov.h"
#include "include/HYPRE.h"
#include "include/HYPRE_parcsr_ls.h"

#include <vector>
#include <numeric>

using namespace espreso;

HypreData::HypreData(MPI_Comm &comm, esint nrows)
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

void HypreData::insertCSR(esint nrows, esint offset, esint *rowPrts, esint *colIndices, double *values, double *rhsValues)
{
	std::vector<esint> ncols, rows;
	ncols.reserve(nrows);
	rows.reserve(nrows);
	for (esint r = 0; r < nrows; r++) {
		ncols.push_back(rowPrts[r + 1] - rowPrts[r]);
		rows.push_back(_roffset + offset + r + 1);
	}
	std::vector<double> x(nrows);

	HYPRE_IJMatrixSetValues(_K, nrows, ncols.data(), rows.data(), colIndices, values);
	HYPRE_IJVectorSetValues(_f, nrows, rows.data(), rhsValues);
	HYPRE_IJVectorSetValues(_x, nrows, rows.data(), x.data());

	_finalized = false;
}

void HypreData::insertIJV(esint nrows, esint offset, esint size, esint *rowIndices, esint *colIndices, double *values, double *rhsValues)
{
	std::vector<esint> ncols, rows;
	ncols.reserve(nrows);
	rows.reserve(nrows);
	for (esint i = 0; i < size; i++) {
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

void HYPRE::solve(const MultigridConfiguration &configuration, HypreData &data, esint nrows, double *solution)
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
	switch (configuration.solver) {
	case MultigridConfiguration::SOLVER::CG:
		HYPRE_ParCSRPCGCreate(data._comm, &solver);

		HYPRE_PCGSetMaxIter(solver, configuration.max_iterations);
		HYPRE_PCGSetTol(solver, configuration.precision);
		HYPRE_PCGSetTwoNorm(solver, 1);
		HYPRE_PCGSetPrintLevel(solver, 0);
		HYPRE_PCGSetLogging(solver, 0);
		break;
	case MultigridConfiguration::SOLVER::GMRES:
	case MultigridConfiguration::SOLVER::FGMRES:
	case MultigridConfiguration::SOLVER::BICGS:
	case MultigridConfiguration::SOLVER::BICGSTAB:
	case MultigridConfiguration::SOLVER::TFQMR:
	case MultigridConfiguration::SOLVER::SYMQMR:
	case MultigridConfiguration::SOLVER::SUPERLU:
	case MultigridConfiguration::SOLVER::SUPERLUX:
	default:
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: not implemented interface to the required solver.";
	}

	HYPRE_Solver preconditioner;
	switch (configuration.preconditioner) {
	case MultigridConfiguration::PRECONDITIONER::DIAGONAL:
	case MultigridConfiguration::PRECONDITIONER::PILUT:
	case MultigridConfiguration::PRECONDITIONER::EUCLID:
	case MultigridConfiguration::PRECONDITIONER::PARASAILS:
	case MultigridConfiguration::PRECONDITIONER::BOOMERAMG:
		HYPRE_BoomerAMGCreate(&preconditioner);
		HYPRE_BoomerAMGSetPrintLevel(preconditioner, 0);
		HYPRE_BoomerAMGSetCoarsenType(preconditioner, 6);
		HYPRE_BoomerAMGSetOldDefault(preconditioner);
		HYPRE_BoomerAMGSetRelaxType(preconditioner, 6);
		HYPRE_BoomerAMGSetNumSweeps(preconditioner, 1);
		HYPRE_BoomerAMGSetTol(preconditioner, 0.0);
		HYPRE_BoomerAMGSetMaxIter(preconditioner, 1);

		HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
		break;
	case MultigridConfiguration::PRECONDITIONER::POLY:
	case MultigridConfiguration::PRECONDITIONER::MLI:
	default:
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: not implemented interface to the required solver.";
	}

	if (run::ecf->environment.print_matrices) {
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::BLUE << "STORE HYPRE SYSTEM";
		{
			std::string prefix = Logging::prepareFile("HYPRE.K");
			HYPRE_IJMatrixPrint(data._K, prefix.c_str());
		}
		{
			std::string prefix = Logging::prepareFile("HYPRE.f");
			HYPRE_IJVectorPrint(data._f, prefix.c_str());
		}
	}

	HYPRE_ParCSRPCGSetup(solver, K, f, x);
	HYPRE_ParCSRPCGSolve(solver, K, f, x);

	esint iterations;
	double norm;
	HYPRE_PCGGetNumIterations(solver, &iterations);
	HYPRE_PCGGetFinalRelativeResidualNorm(solver, &norm);

	ESINFO(CONVERGENCE) << "Final Relative Residual Norm " << norm << " in " << iterations << " iteration.";

	if (run::ecf->environment.print_matrices) {
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::BLUE << "STORE HYPRE SYSTEM SOLUTION";
		std::string prefix = Logging::prepareFile("HYPRE.x");
		HYPRE_IJVectorPrint(data._x, prefix.c_str());
	}

	std::vector<esint> rows(nrows);
	std::iota(rows.begin(), rows.end(), data._roffset + 1);
	HYPRE_IJVectorGetValues(data._x, data._nrows, rows.data(), solution);
}
