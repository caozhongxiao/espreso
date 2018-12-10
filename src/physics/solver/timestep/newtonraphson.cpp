
#include "../../solver/timestep/newtonraphson.h"

#include "../../../config/ecf/physics/physicssolver/nonlinearsolver.h"
#include "../../../basis/logging/logging.h"
#include "../../../linearsolver/linearsolver.h"

#include <cmath>

#include "../../../globals/time.h"
#include "../../dataholder.h"
#include "../../provider/provider.h"
#include "../../solver/loadstep/loadstepsolver.h"

using namespace espreso;

NewtonRaphson::NewtonRaphson(Provider &composer, const NonLinearSolverConfiguration &configuration)
: TimeStepSolver("Newton Raphson", composer), _configuration(configuration)
{

}

void NewtonRaphson::solve(LoadStepSolver &loadStepSolver)
{
	if (!_configuration.check_first_residual && !_configuration.check_second_residual) {
		ESINFO(GLOBAL_ERROR) << "Turn on at least one convergence parameter for NONLINEAR solver.";
	}

	_composer.nextTime();

	Matrices updatedMatrices;
	double &solverPrecision = _composer.linearSolver.precision();
	double solverPrecisionError = 1;

	double temperatureResidual = 10 * _configuration.requested_first_residual;
	double temperatureResidual_first = 0;
	double temperatureResidual_second = 0;

	double heatResidual;
	double heatResidual_first = 0;
	double heatResidual_second = 0;

	double alpha, maxSolutionValue;


	time::iteration = 0;
//	_composer.step.tangentMatrixCorrection = false;
	_composer.solve(loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::B1));
	_composer.parametersChanged();
	_composer.processSolution();
	_composer.storeSubSolution();

//	_composer.step.tangentMatrixCorrection = _configuration.tangent_matrix_correction;
	while (time::iteration++ < _configuration.max_iterations) {
		if (!_configuration.check_second_residual) {
			ESINFO(CONVERGENCE) << "\n >> EQUILIBRIUM ITERATION " << time::iteration + 1 << " IN SUBSTEP "  << time::substep + 1;
		}

		_solution = _composer.instance.primalSolution;
		if (_configuration.method == NonLinearSolverConfiguration::METHOD::NEWTON_RAPHSON) {
			updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::R);
		} else {
			updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::f | Matrices::R);
		}
		if (_configuration.line_search) {
			_f_ext = _composer.instance.f;
		}
		if (_configuration.check_second_residual) {
//			heatResidual_second = _composer.sumSquares(_composer.instance.f, SumRestriction::NON_DIRICHLET, "norm of f not on DIRICHLET");
//			heatResidual_second += _composer.sumSquares(_composer.instance.R, SumRestriction::DIRICHLET, "norm of R on DIRICHLET");
			heatResidual_second = sqrt(heatResidual_second);
			if (heatResidual_second < 1e-3) {
				heatResidual_second = 1e-3;
			}
		}

		_composer.sum(
				_composer.instance.f,
				1, _composer.instance.f,
				-1, _composer.instance.R,
				"f = f - R");

		if (_configuration.check_second_residual) {
			_composer.sum(
					_f_R_BtLambda,
					1, _composer.instance.f,
					-1, _composer.instance.dualSolution,
					"(f - R) - Bt * Lambda");

//			heatResidual_first = sqrt(_composer.sumSquares(_f_R_BtLambda, SumRestriction::NONE, "norm of (f - R) - Bt * Lambda"));
			heatResidual = heatResidual_first / heatResidual_second;

			if (heatResidual < _configuration.requested_second_residual && time::iteration > 1 ) {
				ESINFO(CONVERGENCE) << "    HEAT_CONVERGENCE_VALUE =  " <<  heatResidual_first << "  CRITERION_VALUE = " << heatResidual_second * _configuration.requested_second_residual << " <<< CONVERGED >>>";
				if (_configuration.check_first_residual) {
					if (temperatureResidual < _configuration.requested_first_residual) {
						break;
					}
				} else {
					break;
				}
			} else {
				ESINFO(CONVERGENCE) << " >> EQUILIBRIUM ITERATION " << time::iteration + 1 << " IN SUBSTEP "  << time::substep + 1;
				ESINFO(CONVERGENCE) << "    HEAT_CONVERGENCE_VALUE =  " <<  heatResidual_first << "  CRITERION_VALUE = " << heatResidual_second * _configuration.requested_second_residual;
			}
		}

		if (updatedMatrices & Matrices::K) {
			updatedMatrices |= loadStepSolver.reassembleStructuralMatrices(Matrices::B1c | Matrices::B1duplicity);
		} else {
			updatedMatrices |= loadStepSolver.reassembleStructuralMatrices(Matrices::B1c);
		}
		_composer.addToDirichletInB1(-1, _composer.instance.primalSolution);

		if (_configuration.adaptive_precision) {
			if (time::iteration > 1) {
				solverPrecisionError = temperatureResidual_first / temperatureResidual_second;
				solverPrecision = std::min(_configuration.r_tol * solverPrecisionError, _configuration.c_fact * solverPrecision);
			}
			ESINFO(CONVERGENCE) << "    ADAPTIVE PRECISION = " << solverPrecision << " EPS_ERR = " << solverPrecisionError;
		}

		_composer.solve(updatedMatrices);
		ESINFO(CONVERGENCE) <<  "    LINEAR_SOLVER_OUTPUT: SOLVER = " << "PCG" <<   " N_MAX_ITERATIONS = " << "1" << "  " ;

		if (_configuration.line_search) {
			maxSolutionValue =_composer.maxAbsValue(_composer.instance.primalSolution, "max = |solution|");
			alpha = _composer.lineSearch(_solution, _composer.instance.primalSolution, _f_ext);
			ESINFO(CONVERGENCE) << "    LINE_SEARCH_OUTPUT: " << "PARAMETER = " << alpha << "  MAX_DOF_INCREMENT = " << maxSolutionValue << "  SCALED_MAX_INCREMENT = " << alpha * maxSolutionValue;
		}
		if (_configuration.check_first_residual) {
//			temperatureResidual_first = sqrt(_composer.sumSquares(_composer.instance.primalSolution, SumRestriction::NONE, "|delta U|"));
		}
		_composer.sum(
				_composer.instance.primalSolution,
				1, _composer.instance.primalSolution,
				1, _solution, "U = delta U + U");

		if (_configuration.check_first_residual) {
//			temperatureResidual_second = sqrt(_composer.sumSquares(_composer.instance.primalSolution, SumRestriction::NONE, "|U|"));
			if (temperatureResidual_second < 1e-3) {
				temperatureResidual_second = 1e-3;
			}
			temperatureResidual = temperatureResidual_first / temperatureResidual_second;

			if ( temperatureResidual > _configuration.requested_first_residual){
				ESINFO(CONVERGENCE) << "    TEMPERATURE_CONVERGENCE_VALUE =  " <<  temperatureResidual_first << "  CRITERION_VALUE = " << temperatureResidual_second * _configuration.requested_first_residual ;
			} else {
				ESINFO(CONVERGENCE) << "    TEMPERATURE_CONVERGENCE_VALUE =  " <<  temperatureResidual_first << "  CRITERION_VALUE = " << temperatureResidual_second * _configuration.requested_first_residual <<  " <<< CONVERGED >>>" ;
				if (!_configuration.check_second_residual){
					break;
				}
			}
		}

		_composer.parametersChanged();
		_composer.processSolution();
		_composer.storeSubSolution();
	}

	if (_configuration.check_second_residual) {
		ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << time::iteration;
	} else {
		ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << time::iteration + 1;
	}

	_composer.parametersChanged();
	_composer.processSolution();
	_composer.storeSubSolution();
}


