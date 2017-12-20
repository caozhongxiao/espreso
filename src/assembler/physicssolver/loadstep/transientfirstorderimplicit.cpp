
#include "transientfirstorderimplicit.h"

#include "../assembler.h"
#include "../timestep/timestepsolver.h"
#include "../../step.h"
#include "../../instance.h"
#include "../../physics/physics.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"

#include "../../../basis/logging/logging.h"
#include "../../../config/ecf/physics/physicssolver/transientsolver.h"
#include "../../../config/ecf/environment.h"

#include <iostream>
#include <cmath>

using namespace espreso;

size_t TransientFirstOrderImplicit::loadStep = 0;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(TimeStepSolver &timeStepSolver, const TransientSolverConfiguration &configuration, double duration)
: LoadStepSolver("TRANSIENT", timeStepSolver, duration), _configuration(configuration), _alpha(0), _nTimeStep(_configuration.time_step)
{
	if (configuration.time_step < 1e-7) {
		ESINFO(GLOBAL_ERROR) << "Set time step for TRANSIENT solver greater than 1e-7.";
	}

	U = _assembler.mesh.nodes->appendData({});
	dU = _assembler.mesh.nodes->appendData({});
	V = _assembler.mesh.nodes->appendData({});
	X = _assembler.mesh.nodes->appendData({});
	Y = _assembler.mesh.nodes->appendData({});
	dTK = _assembler.mesh.nodes->appendData({});
}

Matrices TransientFirstOrderImplicit::updateStructuralMatrices(Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::M | Matrices::f | Matrices::R | Matrices::B1 | Matrices::B1c | Matrices::B1duplicity);

	if (_assembler.step.substep) {
		updatedMatrices &= (Matrices::f | Matrices::B1c);
	}

	return reassembleStructuralMatrices(updatedMatrices);
}

Matrices TransientFirstOrderImplicit::reassembleStructuralMatrices(Matrices matrices)
{
	_assembler.updateMatrices(matrices);
	if (matrices & (Matrices::K | Matrices::M)) {
		_assembler.keepK();
		_assembler.sum(
				_assembler.instance.K,
				1 / (_alpha * _assembler.step.timeStep), _assembler.instance.M,
				"K += (1 / alpha * delta T) * M");
	}

	if (matrices & (Matrices::K | Matrices::M | Matrices::f)) {
		_assembler.sum(
				*X->decomposedData,
				1 / (_alpha * _assembler.step.timeStep), *U->decomposedData,
				(1 - _alpha) / _alpha, *V->decomposedData,
				"x = (1 / alpha * delta T) * U + (1 - alpha) / alpha * V");

		_assembler.multiply(
				*Y->decomposedData,
				_assembler.instance.M, *X->decomposedData,
				"y = M * x");

		_assembler.sum(_assembler.instance.f,
				1, _assembler.instance.f,
				1, *Y->decomposedData,
				"f += y");
	}

	return matrices;
}

void TransientFirstOrderImplicit::initLoadStep()
{
	LoadStepSolver::initLoadStep();

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

	if (loadStep + 1 != _assembler.step.step) {
		for (size_t i = 0; i < V->decomposedData->size(); i++) {
			std::fill((*V->decomposedData)[i].begin(), (*V->decomposedData)[i].end(), 0);
		}
	}
	loadStep = _assembler.step.step;
	(*U->decomposedData) = _assembler.instance.primalSolution;
}

void TransientFirstOrderImplicit::runNextTimeStep()
{
	double last = _assembler.step.currentTime;
	_assembler.step.currentTime += _nTimeStep;
	if (_assembler.step.currentTime + _precision >= _startTime + _duration) {
		_assembler.step.currentTime = _startTime + _duration;
	}
	_assembler.step.timeStep = _assembler.step.currentTime - last;
	processTimeStep();
}

void TransientFirstOrderImplicit::processTimeStep()
{
	_assembler.step.internalForceReduction = 1;
	_assembler.step.timeIntegrationConstantK = 1;
	_assembler.step.timeIntegrationConstantM = 1 / (_alpha * _assembler.step.timeStep);

	_timeStepSolver.solve(*this);

	_assembler.sum(
			*dU->decomposedData,
			1, _assembler.instance.primalSolution,
			-1, *U->decomposedData,
			"delta U = U_i - U_i_1");

	_nTimeStep = _assembler.step.timeStep;

	if (_configuration.auto_time_stepping.allowed && _assembler.step.currentTime < _startTime + _duration) {
		double resFreq, oscilationLimit;

		double norm =
				sqrt(_assembler.sumSquares(*dU->decomposedData, SumRestriction::NONE, "|dU|")) /
				sqrt(_assembler.sumSquares(*U->decomposedData, SumRestriction::NONE, "|U|"));

		if (norm < 1e-5) {
			_nTimeStep = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * _assembler.step.timeStep);
		} else {
			_assembler.multiply(
					*dTK->decomposedData,
					_assembler.instance.origK,
					*dU->decomposedData,
					"dTK = delta T * K");
			double TKT = _assembler.multiply(*dTK->decomposedData, *dU->decomposedData, "res. freq. = dTK * delta T");

			_assembler.multiply(
					*dTK->decomposedData,
					_assembler.instance.M,
					*dU->decomposedData,
					"dTM = delta T * M");
			double TMT = _assembler.multiply(*dTK->decomposedData, *dU->decomposedData, "res. freq. = dTM * delta T");


			double gTKT, gTMT;
			MPI_Allreduce(&TKT, &gTKT, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
			MPI_Allreduce(&TMT, &gTMT, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);

			resFreq = gTKT / gTMT;

			oscilationLimit = _assembler.step.timeStep * resFreq;

			double t1 = _configuration.auto_time_stepping.oscilation_limit / resFreq;

			if (_assembler.step.timeStep != t1) {
				if (_assembler.step.timeStep < t1) {
					if (_configuration.auto_time_stepping.IDFactor * _assembler.step.timeStep < t1) {
						_nTimeStep = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * _assembler.step.timeStep);
					}
				} else {
					if (_assembler.step.timeStep / _configuration.auto_time_stepping.IDFactor > t1) {
						_nTimeStep = std::max(_configuration.auto_time_stepping.min_time_step, _assembler.step.timeStep / _configuration.auto_time_stepping.IDFactor);
					}
				}
			}

			ESINFO(CONVERGENCE) << "AUTOMATIC TIME STEPPING INFO: RESPONSE EIGENVALUE(" << resFreq << "), OSCILLATION LIMIT(" << oscilationLimit << ")";
		}

		if (std::fabs(_assembler.step.timeStep - _nTimeStep) / _assembler.step.timeStep < _precision) {
			ESINFO(CONVERGENCE) << "TIME STEP UNCHANGED (" << _nTimeStep << ")";
		} else {
			ESINFO(CONVERGENCE) << "NEW TIME STEP " << (_assembler.step.timeStep < _nTimeStep ? "INCREASED " : "DECREASED ") << "TO VALUE: " << _nTimeStep;
		}
	}

	if (_assembler.step.timeStep - _precision < _nTimeStep) {
		_assembler.sum(
				*V->decomposedData,
				1 / (_alpha * _assembler.step.timeStep), *dU->decomposedData,
				- (1 - _alpha) / _alpha, *V->decomposedData,
				"V = (1 / alpha * delta T) * delta U - (1 - alpha) / alpha * V");

		*U->decomposedData = _assembler.instance.primalSolution;
		_assembler.processSolution();
		_assembler.storeSolution();
	} else {
		_assembler.step.currentTime -= _assembler.step.timeStep;
		--_assembler.step.substep;
	}
}


