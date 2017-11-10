
#include "transientfirstorderimplicit.h"

#include "../assembler.h"
#include "../timestep/timestepsolver.h"
#include "../../step.h"
#include "../../instance.h"
#include "../../solution.h"
#include "../../physics/physics.h"
#include "../../../old/mesh/structures/elementtypes.h"
#include "../../../basis/logging/logging.h"
#include "../../../config/ecf/physics/physicssolver/transientsolver.h"
#include "../../../config/ecf/environment.h"

#include <iostream>
#include <cmath>

using namespace espreso;

size_t TransientFirstOrderImplicit::loadStep = 0;
std::vector<Solution*> TransientFirstOrderImplicit::solutions;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(TimeStepSolver &timeStepSolver, const TransientSolverConfiguration &configuration, double duration)
: LoadStepSolver("TRANSIENT", timeStepSolver, duration), _configuration(configuration), _alpha(0), _nTimeStep(_configuration.time_step)
{
	if (configuration.time_step < 1e-7) {
		ESINFO(GLOBAL_ERROR) << "Set time step for TRANSIENT solver greater than 1e-7.";
	}
}

Matrices TransientFirstOrderImplicit::updateStructuralMatrices(Step &step, Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::M | Matrices::f | Matrices::R | Matrices::B1 | Matrices::B1c | Matrices::B1duplicity);

	if (step.substep && !_timeDependent && !_tempDependent) {
		updatedMatrices &= (Matrices::f | Matrices::B1c);
	}

	return reassembleStructuralMatrices(step, updatedMatrices);
}

Matrices TransientFirstOrderImplicit::reassembleStructuralMatrices(Step &step, Matrices matrices)
{
	_assembler.updateMatrices(step, matrices);
	if (matrices & (Matrices::K | Matrices::M)) {
		_assembler.keepK(step);
		_assembler.sum(
				_assembler.instance.K,
				1 / (_alpha * step.timeStep), _assembler.instance.M,
				"K += (1 / alpha * delta T) * M");
	}

	if (matrices & (Matrices::K | Matrices::M | Matrices::f)) {
		_assembler.sum(
				solutions[SolutionIndex::X]->data,
				1 / (_alpha * step.timeStep), solutions[SolutionIndex::U]->data,
				(1 - _alpha) / _alpha, solutions[SolutionIndex::V]->data,
				"x = (1 / alpha * delta T) * U + (1 - alpha) / alpha * V");

		_assembler.multiply(
				solutions[SolutionIndex::Y]->data,
				_assembler.instance.M, solutions[SolutionIndex::X]->data,
				"y = M * x");

		_assembler.sum(_assembler.instance.f,
				1, _assembler.instance.f,
				1, solutions[SolutionIndex::Y]->data,
				"f += y");
	}

	return matrices;
}

void TransientFirstOrderImplicit::initLoadStep(Step &step)
{
	LoadStepSolver::initLoadStep(step);

	_assembler.setEmptyRegularizationCallback();
	_assembler.setRegularizationFromOrigKCallback();
	_assembler.setB0Callback();

	switch (_configuration.method) {
	case TransientSolverConfiguration::METHOD::CRANK_NICOLSON:
		_alpha = 0.5;
		break;
	case TransientSolverConfiguration::METHOD::GALERKIN:
		_alpha = 2 / 3;
		break;
	case TransientSolverConfiguration::METHOD::BACKWARD_DIFF:
		_alpha = 1;
		break;
	case TransientSolverConfiguration::METHOD::USER:
		_alpha = _configuration.alpha;
		if (_alpha <= 0 || _alpha > 1) {
			ESINFO(GLOBAL_ERROR) << "Alpha has to be from interval (0, 1>.";
		}
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not supported first order implicit solver method.";
	}

	if (!solutions.size()) {
		solutions.push_back(_assembler.addSolution("trans_U" , ElementType::NODES));
		solutions.push_back(_assembler.addSolution("trans_dU", ElementType::NODES));
		solutions.push_back(_assembler.addSolution("trans_V" , ElementType::NODES));
		solutions.push_back(_assembler.addSolution("trans_X" , ElementType::NODES));
		solutions.push_back(_assembler.addSolution("trans_Y" , ElementType::NODES));
		solutions.push_back(_assembler.addSolution("dTK"     , ElementType::NODES));
	}
	if (loadStep + 1 != step.step) {
		solutions[SolutionIndex::V]->fill(0);
	}
	loadStep = step.step;
	solutions[SolutionIndex::U]->data = _assembler.instance.primalSolution;
}

void TransientFirstOrderImplicit::runNextTimeStep(Step &step)
{
	double last = step.currentTime;
	step.currentTime += _nTimeStep;
	if (step.currentTime + _precision >= _startTime + _duration) {
		step.currentTime = _startTime + _duration;
	}
	step.timeStep = step.currentTime - last;
	processTimeStep(step);
}

void TransientFirstOrderImplicit::processTimeStep(Step &step)
{
	step.internalForceReduction = 1;
	step.timeIntegrationConstantK = 1;
	step.timeIntegrationConstantM = 1 / (_alpha * step.timeStep);

	_timeStepSolver.solve(step, *this);

	_assembler.sum(
			solutions[SolutionIndex::dU]->data,
			1, _assembler.instance.primalSolution,
			-1, solutions[SolutionIndex::U]->data,
			"delta U = U_i - U_i_1");

	_nTimeStep = step.timeStep;

	if (_configuration.auto_time_stepping.allowed && step.currentTime < _startTime + _duration) {
		double resFreq, oscilationLimit;

		double norm =
				sqrt(_assembler.sumSquares(step, solutions[SolutionIndex::dU]->data, SumOperation::AVERAGE, SumRestriction::NONE, "|dU|")) /
				sqrt(_assembler.sumSquares(step, solutions[SolutionIndex::U]->data, SumOperation::AVERAGE, SumRestriction::NONE, "|U|"));

		if (norm < 1e-5) {
			_nTimeStep = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * step.timeStep);
		} else {
			_assembler.multiply(
					solutions[SolutionIndex::dTK]->data,
					_assembler.instance.origK,
					solutions[SolutionIndex::dU]->data,
					"dTK = delta T * K");
			double TKT = _assembler.multiply(solutions[SolutionIndex::dTK]->data, solutions[SolutionIndex::dU]->data, "res. freq. = dTK * delta T");

			_assembler.multiply(
					solutions[SolutionIndex::dTK]->data,
					_assembler.instance.M,
					solutions[SolutionIndex::dU]->data,
					"dTM = delta T * M");
			double TMT = _assembler.multiply(solutions[SolutionIndex::dTK]->data, solutions[SolutionIndex::dU]->data, "res. freq. = dTM * delta T");


			double gTKT, gTMT;
			MPI_Allreduce(&TKT, &gTKT, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
			MPI_Allreduce(&TMT, &gTMT, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);

			resFreq = gTKT / gTMT;

			oscilationLimit = step.timeStep * resFreq;

			double t1 = _configuration.auto_time_stepping.oscilation_limit / resFreq;

			if (step.timeStep != t1) {
				if (step.timeStep < t1) {
					if (_configuration.auto_time_stepping.IDFactor * step.timeStep < t1) {
						_nTimeStep = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * step.timeStep);
					}
				} else {
					if (step.timeStep / _configuration.auto_time_stepping.IDFactor > t1) {
						_nTimeStep = std::max(_configuration.auto_time_stepping.min_time_step, step.timeStep / _configuration.auto_time_stepping.IDFactor);
					}
				}
			}

			ESINFO(CONVERGENCE) << "AUTOMATIC TIME STEPPING INFO: RESPONSE EIGENVALUE(" << resFreq << "), OSCILLATION LIMIT(" << oscilationLimit << ")";
		}

		if (std::fabs(step.timeStep - _nTimeStep) / step.timeStep < _precision) {
			ESINFO(CONVERGENCE) << "TIME STEP UNCHANGED (" << _nTimeStep << ")";
		} else {
			ESINFO(CONVERGENCE) << "NEW TIME STEP " << (step.timeStep < _nTimeStep ? "INCREASED " : "DECREASED ") << "TO VALUE: " << _nTimeStep;
		}
	}

	if (step.timeStep - _precision < _nTimeStep) {
		_assembler.sum(
				solutions[SolutionIndex::V]->data,
				1 / (_alpha * step.timeStep), solutions[SolutionIndex::dU]->data,
				- (1 - _alpha) / _alpha, solutions[SolutionIndex::V]->data,
				"V = (1 / alpha * delta T) * delta U - (1 - alpha) / alpha * V");

		solutions[SolutionIndex::U]->data = _assembler.instance.primalSolution;
		_assembler.processSolution(step);
		_assembler.storeSolution(step);
	} else {
		step.currentTime -= step.timeStep;
		--step.substep;
	}
}


