
#include "newtonraphson.h"

#include "esinfo/meshinfo.h"
#include "esinfo/timeinfo.h"
#include "esinfo/eslog.hpp"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/assembler.h"
#include "physics/assembler/composer/composer.h"
#include "physics/assembler/controllers/controller.h"
#include "physics/solver/loadstep/loadstepsolver.h"

#include "config/ecf/physics/physicssolver/nonlinearsolver.h"
#include "basis/containers/serializededata.h"

#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/statisticsstore.h"
#include "linearsolver/linearsolver.h"

#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

NewtonRaphson::NewtonRaphson(NewtonRaphson *previous, Assembler &assembler, NonLinearSolverConfiguration &configuration)
: TimeStepSolver(assembler), _configuration(configuration)
{
	if (previous) {
		_solution = previous->_solution;
	} else {
		_solution = info::mesh->nodes->appendData(_assembler.controller()->solution()->dimension, {});
	}
}

bool NewtonRaphson::hasSameMode(const LoadStepConfiguration &configuration) const
{
	return configuration.mode == LoadStepConfiguration::MODE::NONLINEAR;
}

std::string NewtonRaphson::name()
{
	return "NEWTON RAPHSON";
}

void NewtonRaphson::solve(LoadStepSolver &loadStepSolver)
{
	if (!_configuration.check_first_residual && !_configuration.check_second_residual) {
		eslog::globalerror("Turn on at least one convergence parameter for NONLINEAR solver.\n");
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
			eslog::solver("EQUILIBRIUM ITERATION %d IN SUBSTEP %d\n", time::iteration + 1, time::substep + 1);
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
				eslog::solver("HEAT_CONVERGENCE_VALUE = %.3e / *.3e <<< CONVERGED >>>\n", residualNumerator, residualDenominator * _configuration.requested_second_residual);
				if (_configuration.check_first_residual) {
					if (solutionNorm < _configuration.requested_first_residual) {
						break;
					}
				} else {
					break;
				}
			} else {
				eslog::solver("EQUILIBRIUM ITERATION %d IN SUBSTEP %d\n", time::iteration + 1, time::substep + 1);
				eslog::solver("HEAT_CONVERGENCE_VALUE = %.3e / *.3e\n", residualNumerator, residualDenominator * _configuration.requested_second_residual);
			}
		}

		_assembler.setDirichlet(updatedMatrices, _solution->data);

		if (_configuration.adaptive_precision) {
			double solutionPrecisionError = 1;
			if (time::iteration > 1) {
				solutionPrecisionError = solutionNumerator / solutionDenominator;
				solutionPrecision = std::min(_configuration.r_tol * solutionPrecisionError, _configuration.c_fact * solutionPrecision);
			}
			eslog::solver("ADAPTIVE PRECISION = %.3e EPS_ERR = %.3e\n", solutionPrecision, solutionPrecisionError);
		}

		_assembler.solve(updatedMatrices | Matrices::Dirichlet);
//		eslog::solver("LINEAR_SOLVER_OUTPUT: SOLVER = PCG, N_MAX_ITERATIONS = 1"\n);

		if (_configuration.line_search) {
			double maxSolutionValue = _assembler.controller()->solution()->maxabs();
			double alpha = _assembler.composer()->lineSearch(_solution, _assembler.parameters);
			eslog::solver("LINE SEARCH OUTPUT: PARAMETER = %.3e, MAX DOF INCREMENT = %.3e, SCALED MAX INCREMENT = %.3e\n", alpha, maxSolutionValue, alpha * maxSolutionValue);
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
				eslog::solver("TEMPERATURE CONVERGENCE VALUE = %.3e, CRITERION VALUE = %.3e\n", solutionNumerator, solutionDenominator * _configuration.requested_first_residual);
			} else {
				eslog::solver("TEMPERATURE CONVERGENCE VALUE = %.3e, CRITERION VALUE = %.3e <<<CONVERGED>>>\n", solutionNumerator, solutionDenominator * _configuration.requested_first_residual);
				if (!_configuration.check_second_residual){
					break;
				}
			}
		}
	}

	if (_configuration.check_second_residual) {
		eslog::solver(" >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION %d\n", time::iteration);
	} else {
		eslog::solver(" >> SOLUTION CONVERGED AFTER EQUILIBRIUM ITERATION %d\n", time::iteration + 1);
	}
}


