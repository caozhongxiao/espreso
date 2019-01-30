
#include "physics/assembler/dataholder.h"
#include "esinfo/time.h"
#include "esinfo/meshinfo.h"
#include "newtonraphson.h"
#include "physics/solver/loadstep/loadstepsolver.h"

#include "physics/assembler/assembler.h"
#include "physics/assembler/composer/composer.h"
#include "physics/assembler/controllers/controller.h"

#include "config/ecf/physics/physicssolver/nonlinearsolver.h"
#include "basis/logging/logging.h"
#include "basis/containers/serializededata.h"

#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/statisticsstore.h"
#include "linearsolver/linearsolver.h"

using namespace espreso;

NewtonRaphson::NewtonRaphson(Assembler &assembler, NonLinearSolverConfiguration &configuration)
: TimeStepSolver(assembler), _configuration(configuration)
{
	_solution = info::mesh->nodes->appendData(_assembler.controller()->solution()->dimension, {});
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

	double &solutionPrecision = _assembler.solutionPrecision();

	double solutionNorm = 10 * _configuration.requested_first_residual;
	double solutionNumerator = 0;
	double solutionDenominator = 0;

	time::iteration = 0;
	_assembler.parameters.tangentMatrixCorrection = false;
	Matrices updatedMatrices = loadStepSolver.updateStructuralMatrices(Matrices::K | Matrices::M | Matrices::f);
	_assembler.setDirichlet(updatedMatrices);
	_assembler.solve(updatedMatrices | Matrices::Dirichlet);
	_assembler.parametersChanged();

	_assembler.parameters.tangentMatrixCorrection = _configuration.tangent_matrix_correction;
	while (time::iteration++ < _configuration.max_iterations) {
		if (!_configuration.check_second_residual) {
			ESINFO(CONVERGENCE) << "\n >> EQUILIBRIUM ITERATION " << time::iteration + 1 << " IN SUBSTEP "  << time::substep + 1;
		}

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
			if (residualNumerator / residualDenominator < _configuration.requested_second_residual && time::iteration > 1 ) {
				ESINFO(CONVERGENCE) << "    HEAT_CONVERGENCE_VALUE =  " <<  residualNumerator << "  CRITERION_VALUE = " << residualDenominator * _configuration.requested_second_residual << " <<< CONVERGED >>>";
				if (_configuration.check_first_residual) {
					if (solutionNorm < _configuration.requested_first_residual) {
						break;
					}
				} else {
					break;
				}
			} else {
				ESINFO(CONVERGENCE) << " >> EQUILIBRIUM ITERATION " << time::iteration + 1 << " IN SUBSTEP "  << time::substep + 1;
				ESINFO(CONVERGENCE) << "    HEAT_CONVERGENCE_VALUE =  " <<  residualNumerator << "  CRITERION_VALUE = " << residualDenominator * _configuration.requested_second_residual;
			}
		}

		_assembler.setDirichlet(updatedMatrices, _solution->data);

		if (_configuration.adaptive_precision) {
			double solutionPrecisionError = 1;
			if (time::iteration > 1) {
				solutionPrecisionError = solutionNumerator / solutionDenominator;
				solutionPrecision = std::min(_configuration.r_tol * solutionPrecisionError, _configuration.c_fact * solutionPrecision);
			}
			ESINFO(CONVERGENCE) << "    ADAPTIVE PRECISION = " << solutionPrecision << " EPS_ERR = " << solutionPrecisionError;
		}

		_assembler.solve(updatedMatrices | Matrices::Dirichlet);
		ESINFO(CONVERGENCE) <<  "    LINEAR_SOLVER_OUTPUT: SOLVER = " << "PCG" <<   " N_MAX_ITERATIONS = " << "1" << "  " ;

		if (_configuration.line_search) {
			double maxSolutionValue = _assembler.controller()->solution()->maxabs();
			double alpha = _assembler.composer()->lineSearch(_solution, _assembler.parameters);
			ESINFO(CONVERGENCE) << "    LINE_SEARCH_OUTPUT: " << "PARAMETER = " << alpha << "  MAX_DOF_INCREMENT = " << maxSolutionValue << "  SCALED_MAX_INCREMENT = " << alpha * maxSolutionValue;
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

			if (solutionNorm > _configuration.requested_first_residual) {
				ESINFO(CONVERGENCE) << "    TEMPERATURE_CONVERGENCE_VALUE =  " <<  solutionNumerator << "  CRITERION_VALUE = " << solutionDenominator * _configuration.requested_first_residual ;
			} else {
				ESINFO(CONVERGENCE) << "    TEMPERATURE_CONVERGENCE_VALUE =  " <<  solutionNumerator << "  CRITERION_VALUE = " << solutionDenominator * _configuration.requested_first_residual <<  " <<< CONVERGED >>>" ;
				if (!_configuration.check_second_residual){
					break;
				}
			}
		}
	}

	if (_configuration.check_second_residual) {
		ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << time::iteration;
	} else {
		ESINFO(CONVERGENCE) <<  " >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION " << time::iteration + 1;
	}
}


