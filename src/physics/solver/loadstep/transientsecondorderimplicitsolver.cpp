
#include "transientsecondorderimplicitsolver.h"

#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/timeinfo.h"
#include "esinfo/eslog.hpp"

#include "basis/evaluator/evaluator.h"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/assembler.h"
#include "physics/assembler/composer/composer.h"
#include "physics/assembler/controllers/controller.h"
#include "physics/solver/timestep/timestepsolver.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"

#include "config/ecf/physics/physicssolver/transientsecondorderimplicit.h"
#include "config/ecf/physics/physicssolver/loadstep.h"

#include <cmath>

using namespace espreso;

TransientSecondOrderImplicit::TransientSecondOrderImplicit(TransientSecondOrderImplicit *previous, Assembler &assembler, TimeStepSolver &timeStepSolver, TransientSecondOrderImplicitSolverConfiguration &configuration, double duration)
: LoadStepSolver(assembler, timeStepSolver, duration),
  _configuration(configuration),
  _alpha(_configuration.alpha), _delta(_configuration.delta),
  _massDamping(0), _stiffnessDamping(0),
  _nTimeShift(_configuration.time_step)
{
	if (configuration.time_step < 1e-7) {
		eslog::globalerror("Set time step for TRANSIENT solver greater than 1e-7.\n");
	}

	if (previous) {
		U = previous->U;
		dU = previous->dU;
		V = previous->V;
		W = previous->W;
		X = previous->X;
		Y = previous->Y;
		Z = previous->Z;
		dTK = previous->dTK;
		dTM = previous->dTM;
	} else {
		int dimension = _assembler.controller()->solution()->dimension;
		U = info::mesh->nodes->appendData(dimension, {});
		dU = info::mesh->nodes->appendData(dimension, {});
		V = info::mesh->nodes->appendData(dimension, {});
		W = info::mesh->nodes->appendData(dimension, {});
		Z = info::mesh->nodes->appendData(dimension, {});
		X = info::mesh->nodes->appendData(dimension, {});
		Y = info::mesh->nodes->appendData(dimension, {});
		dTK = info::mesh->nodes->appendData(dimension, {});
		dTM = info::mesh->nodes->appendData(dimension, {});
	}

	_alpha += _configuration.numerical_damping;
	_delta *= (1 + _configuration.numerical_damping) * (1 + _configuration.numerical_damping);
	updateConstants();
	updateDamping();
}

void TransientSecondOrderImplicit::updateConstants()
{
	_newmarkConsts[0] = 1. / (_delta * _nTimeShift * _nTimeShift);
	_newmarkConsts[1] = _alpha / (_delta * _nTimeShift);
	_newmarkConsts[2] = 1. / (_delta * _nTimeShift);
	_newmarkConsts[3] = 1. / (2 * _delta) - 1;
	_newmarkConsts[4] = _alpha / _delta - 1;
	_newmarkConsts[5] = _nTimeShift / 2 * (_alpha / _delta - 2);
	_newmarkConsts[6] = _nTimeShift * (1 - _alpha);
	_newmarkConsts[7] = _nTimeShift * _alpha;
}

void TransientSecondOrderImplicit::updateDamping()
{
	switch (_configuration.damping) {
	case TransientSecondOrderImplicitSolverConfiguration::DAMPING::NONE:
		break;
	case TransientSecondOrderImplicitSolverConfiguration::DAMPING::DIRECT:
		_configuration.direct_damping.stiffness.evaluator->evalVector(1, 0, NULL, NULL, time::current, &_stiffnessDamping);
		_configuration.direct_damping.mass.evaluator->evalVector(1, 0, NULL, NULL, time::current, &_massDamping);
		break;
	case TransientSecondOrderImplicitSolverConfiguration::DAMPING::DAMPING_RATIO: {
		double ratio, frequency;
		_configuration.ratio_damping.ratio.evaluator->evalVector(1, 0, NULL, NULL, time::current, &ratio);
		_configuration.ratio_damping.frequency.evaluator->evalVector(1, 0, NULL, NULL, time::current, &frequency);
		_stiffnessDamping = 2 * ratio * 2 * M_PI * frequency;
		_massDamping = 2 * ratio / (2 * M_PI * frequency);
	} break;
	}
}

bool TransientSecondOrderImplicit::hasSameType(const LoadStepSolverConfiguration &configuration) const
{
	return configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT;
}

std::string TransientSecondOrderImplicit::name()
{
	return "TRANSIENT SECOND ORDER IMPLICIT";
}

Matrices TransientSecondOrderImplicit::updateStructuralMatrices(Matrices matrices)
{
	_assembler.assemble(matrices);

	switch (_configuration.damping) {
	case TransientSecondOrderImplicitSolverConfiguration::DAMPING::NONE:
		if (matrices & (Matrices::K | Matrices::M)) {
			_assembler.composer()->KplusAlfaM(_newmarkConsts[0] + _newmarkConsts[1] * _massDamping);
		}
		if (matrices & (Matrices::K | Matrices::M | Matrices::f)) {
			_assembler.composer()->sum(X, _newmarkConsts[0], U, _newmarkConsts[2], V);
			_assembler.composer()->sum(X, 1, X, _newmarkConsts[3], W);
			_assembler.composer()->applyM(Y, X);
			_assembler.composer()->enrichRHS(1, Y);
		}
		break;
	case TransientSecondOrderImplicitSolverConfiguration::DAMPING::DIRECT:
	case TransientSecondOrderImplicitSolverConfiguration::DAMPING::DAMPING_RATIO:
		if (matrices & (Matrices::K | Matrices::M)) {
			_assembler.composer()->alfaKplusBetaM(1 + _newmarkConsts[1] * _stiffnessDamping, _newmarkConsts[0] + _newmarkConsts[1] * _massDamping);
		}
		if (matrices & (Matrices::K | Matrices::M | Matrices::f)) {
			_assembler.composer()->sum(X, _newmarkConsts[0], U, _newmarkConsts[2], V);
			_assembler.composer()->sum(X, 1, X, _newmarkConsts[3], W);
			_assembler.composer()->sum(X, 1, X, _massDamping * _newmarkConsts[1], U);
			_assembler.composer()->sum(X, 1, X, _massDamping * _newmarkConsts[4], V);
			_assembler.composer()->sum(X, 1, X, _massDamping * _newmarkConsts[5], W);
			_assembler.composer()->applyM(Y, X);
			_assembler.composer()->enrichRHS(1, Y);

			_assembler.composer()->sum(X, _stiffnessDamping * _newmarkConsts[1], U, _stiffnessDamping * _newmarkConsts[4], V);
			_assembler.composer()->sum(X, 1, X, _stiffnessDamping * _newmarkConsts[5], W);
			_assembler.composer()->applyOriginalK(Y, X);
			_assembler.composer()->enrichRHS(1, Y);
		}
		break;
	}

	_assembler.store("transient", matrices);
	return matrices;
}

void TransientSecondOrderImplicit::initLoadStep()
{
	LoadStepSolver::initLoadStep();

	U->data = _assembler.controller()->solution()->data;
}

void TransientSecondOrderImplicit::runNextTimeStep()
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

void TransientSecondOrderImplicit::processTimeStep()
{
	switch (_configuration.damping) {
	case TransientSecondOrderImplicitSolverConfiguration::DAMPING::NONE:
		break;
	case TransientSecondOrderImplicitSolverConfiguration::DAMPING::DIRECT:
		if (_configuration.direct_damping.stiffness.evaluator->isTimeDependent() || _configuration.direct_damping.mass.evaluator->isTimeDependent()) {
			updateDamping();
		}
		break;
	case TransientSecondOrderImplicitSolverConfiguration::DAMPING::DAMPING_RATIO:
		if (_configuration.ratio_damping.ratio.evaluator->isTimeDependent() || _configuration.ratio_damping.frequency.evaluator->isTimeDependent()) {
			updateDamping();
		}
		break;
	}

	_assembler.parameters.internalForceReduction = 1;
	_assembler.parameters.timeIntegrationConstantK = 1 + _newmarkConsts[1] * _stiffnessDamping;
	_assembler.parameters.timeIntegrationConstantM = _newmarkConsts[0] + _newmarkConsts[1] * _massDamping;

	_timeStepSolver.solve(*this);

	_assembler.composer()->sum(dU, 1, _assembler.controller()->solution(), -1, U);
	_nTimeShift = time::shift;

	bool changeConstants = false;

	if (_configuration.auto_time_stepping.allowed && time::current < _startTime + _duration) {
		if (false && dU->norm() / U->norm() < 1e-5) {
			_nTimeShift = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * time::shift);
		} else {
			_assembler.composer()->applyOriginalK(dTK, dU);
			_assembler.composer()->applyM(dTM, dU);

			double resFreq = std::sqrt(std::fabs(_assembler.composer()->mult(dU, dTK) / (4 * M_PI * M_PI * _assembler.composer()->mult(dU, dTM))));
			double resPeriod = 1 / resFreq;
			double t2 = resPeriod / _configuration.auto_time_stepping.points_per_period;

			if (time::shift != t2) {
				if (time::shift < t2) {
					if (_configuration.auto_time_stepping.IDFactor * time::shift < t2) {
						_nTimeShift = std::min(_configuration.auto_time_stepping.max_time_step, _configuration.auto_time_stepping.IDFactor * time::shift);
					}
				} else {
					if (time::shift / _configuration.auto_time_stepping.IDFactor > t2) {
						_nTimeShift = std::max(_configuration.auto_time_stepping.min_time_step, time::shift / _configuration.auto_time_stepping.IDFactor);
					}
				}
			}

			eslog::solver("AUTOMATIC TIME STEPPING INFO: RESPONSE FREQUENCY(%e), POINTS PER PERIOD(%e)\n", resFreq, resPeriod / time::shift);
		}

		if (std::fabs(time::shift - _nTimeShift) / time::shift < _precision) {
			eslog::solver("TIME STEP UNCHANGED (%e)\n", _nTimeShift);
		} else {
			if (time::shift - _precision < _nTimeShift) {
				eslog::solver("INCREASE TIME STEP (%e)\n", _nTimeShift);
			} else {
				eslog::solver("DECREASE TIME STEP (%e)\n", _nTimeShift);
			}
			changeConstants = true;
		}
	}

	if (time::shift - _precision < _nTimeShift) {
		_assembler.composer()->sum(Z, _newmarkConsts[0], dU, -_newmarkConsts[2], V);
		_assembler.composer()->sum(Z, 1, Z, -_newmarkConsts[3], W);
		_assembler.composer()->sum(V, 1, V, _newmarkConsts[6], W);
		_assembler.composer()->sum(V, 1, V, _newmarkConsts[7], Z);
		std::swap(W, Z);

		U->data = _assembler.controller()->solution()->data;
		_assembler.postProcess();
	} else {
		_assembler.controller()->solution()->data = U->data;
		time::current -= time::shift;
		--time::substep;
	}

	if (changeConstants) {
		updateConstants();
	}
}


