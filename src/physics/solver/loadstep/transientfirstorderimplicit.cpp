
#include "transientfirstorderimplicit.h"
#include "../timestep/timestepsolver.h"

#include "../../dataholder.h"
#include "../../assembler/assembler.h"

#include "../../../globals/run.h"
#include "../../../globals/time.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/nodestore.h"

#include "../../../basis/logging/logging.h"
#include "../../../config/ecf/physics/physicssolver/transientsolver.h"

using namespace espreso;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(Assembler &assembler, TimeStepSolver &timeStepSolver, TransientSolverConfiguration &configuration, double duration)
: LoadStepSolver(assembler, timeStepSolver, duration), _configuration(configuration), _alpha(0), _nTimeStep(_configuration.time_step)
{
	if (configuration.time_step < 1e-7) {
		ESINFO(GLOBAL_ERROR) << "Set time step for TRANSIENT solver greater than 1e-7.";
	}

	U = run::mesh->nodes->appendData(1, {});
	dU = run::mesh->nodes->appendData(1, {});
	V = run::mesh->nodes->appendData(1, {});
	X = run::mesh->nodes->appendData(1, {});
	Y = run::mesh->nodes->appendData(1, {});
	dTK = run::mesh->nodes->appendData(1, {});
}

std::string TransientFirstOrderImplicit::name()
{
	return "TRANSIENT";
}

Matrices TransientFirstOrderImplicit::updateStructuralMatrices(Matrices matrices)
{
	Matrices updatedMatrices = matrices & (Matrices::K | Matrices::M | Matrices::f | Matrices::R | Matrices::Dirichlet);

	if (time::substep) {
		updatedMatrices &= (Matrices::f | Matrices::Dirichlet);
	}

	_assembler.assemble(updatedMatrices);
	if (matrices & (Matrices::K | Matrices::M)) {
//		_assembler.keepK();
//		_composer.sum(
//				_composer.instance.K,
//				1 / (_alpha * time::shift), _composer.instance.M,
//				"K += (1 / alpha * delta T) * M");
	}

//	_assembler.updateGluingMatrices(updatedMatrices);

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

//	(*U->decomposedData) = _composer.instance.primalSolution;
}

void TransientFirstOrderImplicit::runNextTimeStep()
{
	double last = time::current;
	time::current += _nTimeStep;
	if (time::current + _precision >= _startTime + _duration) {
		time::current = _startTime + _duration;
	}
	time::shift = time::current - last;
	processTimeStep();
}

void TransientFirstOrderImplicit::processTimeStep()
{
	_assembler.parameters.internalForceReduction = 1;
	_assembler.parameters.timeIntegrationConstantK = 1;
	_assembler.parameters.timeIntegrationConstantM = 1 / (_alpha * time::shift);

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


