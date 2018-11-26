
#include "../../solver/loadstep/transientfirstorderimplicit.h"

#include "../../step.h"
#include "../../instance.h"
#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"

#include "../../../basis/logging/logging.h"
#include "../../../config/ecf/physics/physicssolver/transientsolver.h"
#include "../../../config/ecf/environment.h"

#include <iostream>
#include <cmath>

#include "../../assembler/physics.h"
#include "../../provider/provider.h"
#include "../../solver/timestep/timestepsolver.h"

using namespace espreso;

size_t TransientFirstOrderImplicit::loadStep = 0;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(TimeStepSolver &timeStepSolver, const TransientSolverConfiguration &configuration, double duration)
: LoadStepSolver("TRANSIENT", timeStepSolver, duration), _configuration(configuration), _alpha(0), _nTimeStep(_configuration.time_step)
{
	if (configuration.time_step < 1e-7) {
		ESINFO(GLOBAL_ERROR) << "Set time step for TRANSIENT solver greater than 1e-7.";
	}

	U = _composer.mesh.nodes->appendData(1, {});
	dU = _composer.mesh.nodes->appendData(1, {});
	V = _composer.mesh.nodes->appendData(1, {});
	X = _composer.mesh.nodes->appendData(1, {});
	Y = _composer.mesh.nodes->appendData(1, {});
	dTK = _composer.mesh.nodes->appendData(1, {});
}

Matrices TransientFirstOrderImplicit::updateStructuralMatrices(Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::M | Matrices::f | Matrices::R | Matrices::B1);

//	if (_composer.step.substep) {
//		updatedMatrices &= (Matrices::f | Matrices::B1c);
//	}

	return reassembleStructuralMatrices(updatedMatrices);
}

Matrices TransientFirstOrderImplicit::reassembleStructuralMatrices(Matrices matrices)
{
	_composer.updateStructuralMatrices(matrices);
	if (matrices & (Matrices::K | Matrices::M)) {
		_composer.keepK();
		_composer.sum(
				_composer.instance.K,
				1 / (_alpha * _composer.step.timeStep), _composer.instance.M,
				"K += (1 / alpha * delta T) * M");
	}

	_composer.updateGluingMatrices(matrices);

	if (matrices & (Matrices::K | Matrices::M | Matrices::f)) {
//		_composer.sum(
//				*X->decomposedData,
//				1 / (_alpha * _composer.step.timeStep), *U->decomposedData,
//				(1 - _alpha) / _alpha, *V->decomposedData,
//				"x = (1 / alpha * delta T) * U + (1 - alpha) / alpha * V");
//
//		_composer.multiply(
//				*Y->decomposedData,
//				_composer.instance.M, *X->decomposedData,
//				"y = M * x");
//
//		_composer.sum(_composer.instance.f,
//				1, _composer.instance.f,
//				1, *Y->decomposedData,
//				"f += y");
	}

	return matrices;
}

void TransientFirstOrderImplicit::initLoadStep()
{
	LoadStepSolver::initLoadStep();

	_composer.setEmptyRegularizationCallback();
	_composer.setRegularizationFromOrigKCallback();
	_composer.setB0Callback();

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

	if (loadStep + 1 != _composer.step.step) {
//		for (size_t i = 0; i < V->decomposedData->size(); i++) {
//			std::fill((*V->decomposedData)[i].begin(), (*V->decomposedData)[i].end(), 0);
//		}
	}
	loadStep = _composer.step.step;
//	(*U->decomposedData) = _composer.instance.primalSolution;
}

void TransientFirstOrderImplicit::runNextTimeStep()
{
	double last = _composer.step.currentTime;
	_composer.step.currentTime += _nTimeStep;
	if (_composer.step.currentTime + _precision >= _startTime + _duration) {
		_composer.step.currentTime = _startTime + _duration;
	}
	_composer.step.timeStep = _composer.step.currentTime - last;
	processTimeStep();
}

void TransientFirstOrderImplicit::processTimeStep()
{
	_composer.step.internalForceReduction = 1;
	_composer.step.timeIntegrationConstantK = 1;
	_composer.step.timeIntegrationConstantM = 1 / (_alpha * _composer.step.timeStep);

	_timeStepSolver.solve(*this);

//	_composer.sum(
//			*dU->decomposedData,
//			1, _composer.instance.primalSolution,
//			-1, *U->decomposedData,
//			"delta U = U_i - U_i_1");
//
//	_nTimeStep = _composer.step.timeStep;
//
//	if (_configuration.auto_time_stepping.allowed && _composer.step.currentTime < _startTime + _duration) {
//		double resFreq, oscilationLimit;
//
//		double norm =
//				sqrt(_composer.sumSquares(*dU->decomposedData, SumRestriction::NONE, "|dU|")) /
//				sqrt(_composer.sumSquares(*U->decomposedData, SumRestriction::NONE, "|U|"));
//
//		if (norm < 1e-5) {
//			_nTimeStep = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * _composer.step.timeStep);
//		} else {
//			_composer.multiply(
//					*dTK->decomposedData,
//					_composer.instance.origK,
//					*dU->decomposedData,
//					"dTK = delta T * K");
//			double TKT = _composer.multiply(*dTK->decomposedData, *dU->decomposedData, "res. freq. = dTK * delta T");
//
//			_composer.multiply(
//					*dTK->decomposedData,
//					_composer.instance.M,
//					*dU->decomposedData,
//					"dTM = delta T * M");
//			double TMT = _composer.multiply(*dTK->decomposedData, *dU->decomposedData, "res. freq. = dTM * delta T");
//
//
//			double gTKT, gTMT;
//			MPI_Allreduce(&TKT, &gTKT, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
//			MPI_Allreduce(&TMT, &gTMT, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
//
//			resFreq = gTKT / gTMT;
//
//			oscilationLimit = _composer.step.timeStep * resFreq;
//
//			double t1 = _configuration.auto_time_stepping.oscilation_limit / resFreq;
//
//			if (_composer.step.timeStep != t1) {
//				if (_composer.step.timeStep < t1) {
//					if (_configuration.auto_time_stepping.IDFactor * _composer.step.timeStep < t1) {
//						_nTimeStep = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * _composer.step.timeStep);
//					}
//				} else {
//					if (_composer.step.timeStep / _configuration.auto_time_stepping.IDFactor > t1) {
//						_nTimeStep = std::max(_configuration.auto_time_stepping.min_time_step, _composer.step.timeStep / _configuration.auto_time_stepping.IDFactor);
//					}
//				}
//			}
//
//			ESINFO(CONVERGENCE) << "AUTOMATIC TIME STEPPING INFO: RESPONSE EIGENVALUE(" << resFreq << "), OSCILLATION LIMIT(" << oscilationLimit << ")";
//		}
//
//		if (std::fabs(_composer.step.timeStep - _nTimeStep) / _composer.step.timeStep < _precision) {
//			ESINFO(CONVERGENCE) << "TIME STEP UNCHANGED (" << _nTimeStep << ")";
//		} else {
//			ESINFO(CONVERGENCE) << "NEW TIME STEP " << (_composer.step.timeStep < _nTimeStep ? "INCREASED " : "DECREASED ") << "TO VALUE: " << _nTimeStep;
//		}
//	}
//
//	if (_composer.step.timeStep - _precision < _nTimeStep) {
//		_composer.sum(
//				*V->decomposedData,
//				1 / (_alpha * _composer.step.timeStep), *dU->decomposedData,
//				- (1 - _alpha) / _alpha, *V->decomposedData,
//				"V = (1 / alpha * delta T) * delta U - (1 - alpha) / alpha * V");
//
//		*U->decomposedData = _composer.instance.primalSolution;
//		_composer.processSolution();
//		_composer.storeSolution();
//	} else {
//		_composer.step.currentTime -= _composer.step.timeStep;
//		--_composer.step.substep;
//	}
}


