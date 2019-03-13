
#include <config/ecf/physics/physicssolver/nonlinear.h>
#include "newtonraphsonsolver.h"
#include "esinfo/meshinfo.h"
#include "esinfo/timeinfo.h"
#include "esinfo/eslog.hpp"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/assembler.h"
#include "physics/assembler/composer/composer.h"
#include "physics/assembler/controllers/controller.h"
#include "physics/solver/loadstep/loadstepsolver.h"

#include "basis/containers/serializededata.h"

#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/statisticsstore.h"
#include "linearsolver/linearsolver.h"

#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

const char* NewtonRaphson::statusNOT = "NOT";
const char* NewtonRaphson::statusYES = "YES";

NewtonRaphson::NewtonRaphson(NewtonRaphson *previous, Assembler &assembler, NonLinearSolverConfiguration &configuration)
: TimeStepSolver(assembler), _configuration(configuration),
  _lsAlpha(1), _lsInc(0), _lsSInc(0),
  _firstConvergenceValue(0), _firstCriterionValue(0),
  _secondConvergenceValue(0), _secondCriterionValue(0),
  _status(statusNOT)
{
	if (previous) {
		_solution = previous->_solution;
	} else {
		_solution = info::mesh->nodes->appendData(_assembler.controller()->solution()->dimension, {});
	}
}

bool NewtonRaphson::hasSameMode(const LoadStepSolverConfiguration &configuration) const
{
	return configuration.mode == LoadStepSolverConfiguration::MODE::NONLINEAR;
}

std::string NewtonRaphson::name()
{
	return "NEWTON RAPHSON";
}

void NewtonRaphson::setSolverParams()
{
	eslog::addsolverparam("ITERATION [IT] -- the non-linear iteration number", "  IT", "%4d", time::iteration);
	if (_configuration.line_search) {
		eslog::addsolverparam("LINE SEARCH ALPHA [LS A]", " LS A", "%5.3f", _lsAlpha);
		eslog::addsolverparam("LINE SEARCH MAX DOF INCREMENT [LS INC]", "  LS INC", "%.2e", _lsInc);
		eslog::addsolverparam("LINE SEARCH SCALED MAX INCREMENT [LS SINC]", " LS SINC", "%.2e", _lsSInc);
	}
	if (_configuration.check_first_residual) {
		eslog::addsolverparam("TEMPERATURE CONVERGENCE VALUE [T CNV]", "    T CNV", "%.3e", _firstConvergenceValue);
		eslog::addsolverparam("TEMPERATURE CRITERION VALUE [T CRT]", "    T CRT", "%.3e", _firstCriterionValue);
	}
	if (_configuration.check_second_residual) {
		eslog::addsolverparam("HEAT CONVERGENCE VALUE [T CNV]", "    H CNV", "%.3e", _secondConvergenceValue);
		eslog::addsolverparam("HEAT CRITERION VALUE [T CRT]", "    H CRT", "%.3e", _secondCriterionValue);
	}
	eslog::addsolverparam("CONVERGENCE STATUS [CNV]", "CNV", "%s", _status);
}

void NewtonRaphson::solve(LoadStepSolver &loadStepSolver)
{
	if (!_configuration.check_first_residual && !_configuration.check_second_residual) {
		eslog::globalerror("Turn on at least one convergence parameter for NONLINEAR solver.\n");
	}

	switch (_configuration.method) {
	case NonLinearSolverConfiguration::METHOD::NEWTON_RAPHSON:
		eslog::solver("  > ----------------------- NEWTON  RAPHSON ----------------------- <\n");
		break;
	case NonLinearSolverConfiguration::METHOD::MODIFIED_NEWTON_RAPHSON:
		eslog::solver("  > ------------------- MODIFIED NEWTON RAPHSON ------------------- <\n");
		break;
	}

	double &solutionPrecision = _assembler.solutionPrecision();

	double solutionNorm = 10 * _configuration.requested_first_residual;
	double solutionNumerator = 0;
	double solutionDenominator = 0;

	_lsAlpha = 1;
	_lsInc = 0;
	_lsSInc = 0;
	_firstConvergenceValue = 0;
	_firstCriterionValue = 0;
	_secondConvergenceValue = 0;
	_secondCriterionValue = 0;
	_status = statusNOT;
	time::iteration = 0;
	_assembler.parameters.tangentMatrixCorrection = false;
	Matrices updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f);
	_assembler.setDirichlet(updatedMatrices);

	eslog::solver("  > INITIAL STEP              REASSEMBLED MATRICES :: %c, %c, %c, %c, %c <\n",
			(updatedMatrices & Matrices::K) ? 'K' : ' ',
			(updatedMatrices & Matrices::M) ? 'M' : ' ',
			(updatedMatrices & Matrices::C) ? 'C' : ' ',
			(updatedMatrices & Matrices::R) ? 'R' : ' ',
			(updatedMatrices & Matrices::f) ? 'f' : ' ');

	_assembler.solve(updatedMatrices | Matrices::Dirichlet);
	_assembler.parametersChanged();

	_assembler.parameters.tangentMatrixCorrection = _configuration.tangent_matrix_correction;
	while (time::iteration++ < _configuration.max_iterations) {
		_solution->data = _assembler.controller()->solution()->data;
		if (_configuration.method == NonLinearSolverConfiguration::METHOD::NEWTON_RAPHSON) {
			updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::R);
		} else {
			updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::f | Matrices::R);
		}
		_assembler.composer()->keepRHS();
		_assembler.composer()->RHSMinusR();

		if (_configuration.check_second_residual) {
			double residualNumerator = _assembler.composer()->residualNormNumerator();
			double residualDenominator = std::max(_assembler.composer()->residualNormDenominator(), 1e-3);
			_secondConvergenceValue = residualNumerator;
			_secondCriterionValue = residualDenominator * _configuration.requested_second_residual;
			if (residualNumerator / residualDenominator < _configuration.requested_second_residual && time::iteration > 1 ) {
				eslog::solver("  > HEAT NORM, CRITERIA          %.5e, %.5e CONVERGED <\n", residualNumerator, residualDenominator * _configuration.requested_second_residual);
				if (_configuration.check_first_residual) {
					if (solutionNorm < _configuration.requested_first_residual) {
						_status = statusYES;
						break;
					}
				} else {
					_status = statusYES;
					break;
				}
			} else {
				eslog::solver("  > HEAT NORM, CRITERIA          %.5e, %.5e           <\n", residualNumerator, residualDenominator * _configuration.requested_second_residual);
			}
		}

		_assembler.setDirichlet(updatedMatrices, _solution->data);

		if (_configuration.adaptive_precision) {
			double solutionPrecisionError = 1;
			if (time::iteration > 1) {
				solutionPrecisionError = solutionNumerator / solutionDenominator;
				solutionPrecision = std::min(_configuration.r_tol * solutionPrecisionError, _configuration.c_fact * solutionPrecision);
			}
			eslog::solver("  > ADAPTIVE PRECISION, EPS ERROR         %.5e / %.5e <\n", solutionPrecision, solutionPrecisionError);
		}

		eslog::solver("\n  > EQUIL. ITERATION :: %3d   REASSEMBLED MATRICES :: %c, %c, %c, %c, %c <\n",
				time::iteration,
				(updatedMatrices & Matrices::K) ? 'K' : ' ',
				(updatedMatrices & Matrices::M) ? 'M' : ' ',
				(updatedMatrices & Matrices::C) ? 'C' : ' ',
				(updatedMatrices & Matrices::R) ? 'R' : ' ',
				(updatedMatrices & Matrices::f) ? 'f' : ' ');
		_assembler.solve(updatedMatrices | Matrices::Dirichlet);

		if (_configuration.line_search) {
			_lsInc = _assembler.controller()->solution()->maxabs();
			_lsAlpha = _assembler.composer()->lineSearch(_solution, _assembler.parameters);
			_lsSInc = _lsAlpha * _lsInc;
			eslog::solver("  > LINE SEARCH, MAX DOF INCREMENT             %.5f, %.5e <\n", _lsAlpha, _lsInc);
		}

		if (!_configuration.check_first_residual) {
			_assembler.composer()->enrichSolution(1, _solution);
			_assembler.parametersChanged();
		} else {
			solutionNumerator = _assembler.controller()->solution()->norm();
			_assembler.composer()->enrichSolution(1, _solution);
			_assembler.parametersChanged();
			solutionDenominator = std::max(_assembler.controller()->solution()->norm(), 1e-3);
			solutionNorm = solutionNumerator / solutionDenominator;

			_firstConvergenceValue = solutionNumerator;
			_firstCriterionValue = solutionDenominator * _configuration.requested_first_residual;
			if (solutionNorm > _configuration.requested_first_residual) {
				_status = statusNOT;
				eslog::solver("  > TEMPERATURE NORM, CRITERIA  %.5e / %.5e           <\n", solutionNumerator, solutionDenominator * _configuration.requested_first_residual);
			} else {
				_status = statusYES;
				eslog::solver("  > TEMPERATURE NORM, CRITERIA  %.5e / %.5e CONVERGED <\n", solutionNumerator, solutionDenominator * _configuration.requested_first_residual);
				if (!_configuration.check_second_residual) {
					break;
				}
			}
		}
	}

	eslog::solver("  > --------------------------------------------------------------- <\n");
}


