
#include "transientfirstorderimplicitsolver.h"

#include "esinfo/meshinfo.h"
#include "esinfo/timeinfo.h"
#include "esinfo/eslog.hpp"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/assembler.h"
#include "physics/assembler/composer/composer.h"
#include "physics/assembler/controllers/controller.h"
#include "physics/solver/timestep/timestepsolver.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"

#include "config/ecf/physics/physicssolver/transientfirstorderimplicit.h"
#include "config/ecf/physics/physicssolver/loadstep.h"

#include <cmath>

using namespace espreso;

TransientFirstOrderImplicit::TransientFirstOrderImplicit(TransientFirstOrderImplicit *previous, Assembler &assembler, TimeStepSolver &timeStepSolver, TransientFirstOrderImplicitSolverConfiguration &configuration, double duration)
: LoadStepSolver(assembler, timeStepSolver, duration), _configuration(configuration), _alpha(0), _nTimeShift(_configuration.time_step)
{
	if (configuration.time_step < 1e-7) {
		eslog::globalerror("Set time step for TRANSIENT solver greater than 1e-7.\n");
	}

	if (previous) {
		U = previous->U;
		dU = previous->dU;
		V = previous->V;
		X = previous->X;
		Y = previous->Y;
		dTK = previous->dTK;
		dTM = previous->dTM;
	} else {
		U = info::mesh->nodes->appendData(1, {});
		dU = info::mesh->nodes->appendData(1, {});
		V = info::mesh->nodes->appendData(1, {});
		X = info::mesh->nodes->appendData(1, {});
		Y = info::mesh->nodes->appendData(1, {});
		dTK = info::mesh->nodes->appendData(1, {});
		dTM = info::mesh->nodes->appendData(1, {});
	}
}

bool TransientFirstOrderImplicit::hasSameType(const LoadStepSolverConfiguration &configuration) const
{
	return configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT;
}

std::string TransientFirstOrderImplicit::name()
{
	return "TRANSIENT FIRST ORDER IMPLICIT";
}

Matrices TransientFirstOrderImplicit::updateStructuralMatrices(Matrices matrices)
{
	_assembler.assemble(matrices);
	if (matrices & (Matrices::K | Matrices::M)) {
		_assembler.composer()->KplusAlfaM(1 / (_alpha * time::shift));
	}

	if (matrices & (Matrices::K | Matrices::M | Matrices::f)) {
		_assembler.composer()->sum(X, 1 / (_alpha * time::shift), U, (1 - _alpha) / _alpha, V);
		_assembler.composer()->applyM(Y, X);
		_assembler.composer()->enrichRHS(1, Y);
	}

	return matrices;
}

void TransientFirstOrderImplicit::initLoadStep()
{
	LoadStepSolver::initLoadStep();

	switch (_configuration.method) {
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::CRANK_NICOLSON:
		_alpha = 0.5;
		break;
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::GALERKIN:
		_alpha = 2 / 3;
		break;
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::BACKWARD_DIFF:
		_alpha = 1;
		break;
	case TransientFirstOrderImplicitSolverConfiguration::METHOD::USER:
		_alpha = _configuration.alpha;
		if (_alpha <= 0 || _alpha > 1) {
			eslog::globalerror("Alpha has to be from interval (0, 1>.\n");
		}
		break;
	default:
		eslog::globalerror("Not supported first order implicit solver method.\n");
	}

	U->data = _assembler.controller()->solution()->data;
}

void TransientFirstOrderImplicit::runNextTimeStep()
{
	double last = time::current;
	time::current += _nTimeShift;
	if (time::current + _precision >= _startTime + _duration) {
		time::current = _startTime + _duration;
	}
	time::shift = time::current - last;
	_assembler.nextTime();

	processTimeStep();
}

void TransientFirstOrderImplicit::processTimeStep()
{
	_assembler.parameters.internalForceReduction = 1;
	_assembler.parameters.timeIntegrationConstantK = 1;
	_assembler.parameters.timeIntegrationConstantM = 1 / (_alpha * time::shift);

	_timeStepSolver.solve(*this);

	_assembler.composer()->sum(dU, 1, _assembler.controller()->solution(), -1, U);
	_nTimeShift = time::shift;

	if (_configuration.auto_time_stepping.allowed && time::current < _startTime + _duration) {
		if (dU->norm() / U->norm() < 1e-5) {
			_nTimeShift = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * time::shift);
		} else {
			_assembler.composer()->applyOriginalK(dTK, dU);
			_assembler.composer()->applyM(dTM, dU);

			double resFreq = _assembler.composer()->mult(dU, dTK) / _assembler.composer()->mult(dU, dTM);
			double oscilationLimit = time::shift * resFreq;
			double t1 = _configuration.auto_time_stepping.oscilation_limit / resFreq;

			if (time::shift != t1) {
				if (time::shift < t1) {
					if (_configuration.auto_time_stepping.IDFactor * time::shift < t1) {
						_nTimeShift = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * time::shift);
					}
				} else {
					if (time::shift / _configuration.auto_time_stepping.IDFactor > t1) {
						_nTimeShift = std::max(_configuration.auto_time_stepping.min_time_step, time::shift / _configuration.auto_time_stepping.IDFactor);
					}
				}
			}

			eslog::linearsolver("AUTOMATIC TIME STEPPING INFO: RESPONSE EIGENVALUE(%e), OSCILLATION LIMIT(%e)\n", resFreq, oscilationLimit);
		}

		if (std::fabs(time::shift - _nTimeShift) / time::shift < _precision) {
			eslog::linearsolver("TIME STEP UNCHANGED (%5.3f)\n", _nTimeShift);
		} else {
			if (time::shift < _nTimeShift) {
				eslog::linearsolver("INCREASE TIME STEP (%5.3f)\n", _nTimeShift);
			} else {
				eslog::linearsolver("DECREASE TIME STEP (%5.3f)\n", _nTimeShift);
			}
		}
	}

	if (time::shift - _precision < _nTimeShift) {
		_assembler.composer()->sum(V, 1 / (_alpha * time::shift), dU, -(1 - _alpha) / _alpha, V);

		U->data = _assembler.controller()->solution()->data;
		_assembler.postProcess();
	} else {
		_assembler.controller()->solution()->data = U->data;
		time::current -= time::shift;
		--time::substep;
	}
}


