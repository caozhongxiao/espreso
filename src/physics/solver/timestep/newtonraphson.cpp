
#include "newtonraphson.h"
#include "../loadstep/loadstepsolver.h"

#include "../../dataholder.h"
#include "../../assembler/assembler.h"

#include "../../../globals/time.h"
#include "../../../globals/run.h"
#include "../../../config/ecf/physics/physicssolver/nonlinearsolver.h"
#include "../../../basis/logging/logging.h"
#include "../../../basis/containers/serializededata.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"
#include "../../../mesh/store/statisticsstore.h"
#include "../../../linearsolver/linearsolver.h"

using namespace espreso;

NewtonRaphson::NewtonRaphson(Assembler &assembler, LinearSolver &solver, NonLinearSolverConfiguration &configuration)
: TimeStepSolver(assembler, solver), _configuration(configuration)
{
	_solution = run::mesh->nodes->appendData(_assembler.solution()->dimension, {});
	_RHS = run::mesh->nodes->appendData(_assembler.solution()->dimension, {});
}

std::string NewtonRaphson::name()
{
	return "NEWTON RAPHSON";
}

void NewtonRaphson::solve(LoadStepSolver &loadStepSolver)
{
	if (!_configuration.check_first_residual && !_configuration.check_second_residual) {
		ESINFO(GLOBAL_ERROR) << "Turn on at least one convergence parameter for NONLINEAR solver.";
	}

	Matrices updatedMatrices;
	double &solutionPrecision = _assembler.solutionPrecision();

	double temperatureResidual = 10 * _configuration.requested_first_residual;
	double temperatureResidual_first = 0;
	double temperatureResidual_second = 0;

	double heatResidual_first = 0;
	double heatResidual_second = 0;

	std::vector<Statistics> stats(_assembler.solution()->dimension + 1);

	time::iteration = 0;
	_assembler.parameters.tangentMatrixCorrection = false;
	_solver.solve(loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::Dirichlet));

	_assembler.postProcess();

	_assembler.parameters.tangentMatrixCorrection = _configuration.tangent_matrix_correction;
	while (time::iteration++ < _configuration.max_iterations) {
		if (!_configuration.check_second_residual) {
			ESINFO(CONVERGENCE) << "\n >> EQUILIBRIUM ITERATION " << time::iteration + 1 << " IN SUBSTEP "  << time::substep + 1;
		}

		_solution->data = _assembler.solution()->data;
		if (_configuration.method == NonLinearSolverConfiguration::METHOD::NEWTON_RAPHSON) {
			updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f | Matrices::R);
		} else {
			updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::f | Matrices::R);
		}
		if (_configuration.line_search) {
			_RHS = _assembler.RHS();
		}
		if (_configuration.check_second_residual) {
			heatResidual_second = _assembler.residualNorm();
			if (heatResidual_second < 1e-3) {
				heatResidual_second = 1e-3;
			}
		}

		_assembler.RHSMinusR();
		if (_configuration.check_second_residual) {
//			_assembler.sum(
//					_f_R_BtLambda,
//					1, _assembler.instance.f,
//					-1, _assembler.instance.dualSolution,
//					"(f - R) - Bt * Lambda");

//			heatResidual_first = sqrt(_assembler.sumSquares(_f_R_BtLambda, SumRestriction::NONE, "norm of (f - R) - Bt * Lambda"));

			if (heatResidual_first / heatResidual_second < _configuration.requested_second_residual && time::iteration > 1 ) {
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

		updatedMatrices |= loadStepSolver.updateStructuralMatrices(Matrices::Dirichlet);
		_assembler.DirichletMinusRHS();

		if (_configuration.adaptive_precision) {
			double solutionPrecisionError = 1;
			if (time::iteration > 1) {
				solutionPrecisionError = temperatureResidual_first / temperatureResidual_second;
				solutionPrecision = std::min(_configuration.r_tol * solutionPrecisionError, _configuration.c_fact * solutionPrecision);
			}
			ESINFO(CONVERGENCE) << "    ADAPTIVE PRECISION = " << solutionPrecision << " EPS_ERR = " << solutionPrecisionError;
		}

		_solver.solve(updatedMatrices);
		ESINFO(CONVERGENCE) <<  "    LINEAR_SOLVER_OUTPUT: SOLVER = " << "PCG" <<   " N_MAX_ITERATIONS = " << "1" << "  " ;

		if (_configuration.line_search) {
			_assembler.solution()->statistics(run::mesh->allNodes()->nodes->datatarray(), run::mesh->nodes->uniqueTotalSize, stats.data());
			double maxSolutionValue = std::max(std::fabs(stats[0].min), std::fabs(stats[0].max));
//			double alpha = _assembler.lineSearch(_solution, _assembler.instance.primalSolution, _f_ext);
//			ESINFO(CONVERGENCE) << "    LINE_SEARCH_OUTPUT: " << "PARAMETER = " << alpha << "  MAX_DOF_INCREMENT = " << maxSolutionValue << "  SCALED_MAX_INCREMENT = " << alpha * maxSolutionValue;
		}
		if (_configuration.check_first_residual) {
			temperatureResidual_first = _assembler.solution()->norm();
		}
		_assembler.enrichRHS(1, _solution);

		if (_configuration.check_first_residual) {
			temperatureResidual_second = _assembler.solution()->norm();
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

		_assembler.postProcess();
	}

	if (_configuration.check_second_residual) {
		ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << time::iteration;
	} else {
		ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << time::iteration + 1;
	}

	_assembler.postProcess();
}


