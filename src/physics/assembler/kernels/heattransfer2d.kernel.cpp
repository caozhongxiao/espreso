
#include "physics/assembler/dataholder.h"
#include "esinfo/time.h"
#include "esinfo/ecfinfo.h"
#include "heattransfer2d.kernel.h"

#include "physics/assembler/assembler.h"
#include "basis/containers/point.h"
#include "basis/matrices/denseMatrix.h"
#include "basis/evaluator/evaluator.h"

#include "config/ecf/root.h"
#include "mesh/elements/element.h"


using namespace espreso;

// TODO: create file with constants
#define CONST_Stefan_Boltzmann 5.6703e-8

using namespace espreso;

void HeatTransfer2DKernel::assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const
{
	auto d2r = [] (double degree) -> double {
		return M_PI * degree / 180;
	};

	Point p(coordinates[0], coordinates[1], 0);

	double cos, sin;
	switch (mat->coordinate_system.type) {
	case CoordinateSystemConfiguration::TYPE::CARTESIAN:
		cos = std::cos(d2r(mat->coordinate_system.rotation.z.evaluator->evaluate(p, time, temp)));
		sin = std::sin(d2r(mat->coordinate_system.rotation.z.evaluator->evaluate(p, time, temp)));
		break;
	case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: {
		Point origin(
				mat->coordinate_system.center.x.evaluator->evaluate(p, time, temp),
				mat->coordinate_system.center.y.evaluator->evaluate(p, time, temp),
				0);
		double rotation = std::atan2((p.y - origin.y), (p.x - origin.x));
		cos = std::cos(rotation);
		sin = std::sin(rotation);
		break;
	}
	default:
		ESINFO(ERROR) << "Invalid material type (SPHERICAL for 2D).";
	}

	DenseMatrix TCT(2, 2), T(2, 2), C(2, 2), _CD, TCDT;
	T(0, 0) =  cos; T(0, 1) = sin;
	T(1, 0) = -sin; T(1, 1) = cos;

	if (tangentCorrection) {
		_CD.resize(2, 2);
		TCDT.resize(2, 2);
	}

	auto derivation = [&] (const ECFExpression &expression, double h) {
		return (
				expression.evaluator->evaluate(p, time, temp + h) -
				expression.evaluator->evaluate(p, time, temp - h)
				) / (2 * h);
	};

	switch (mat->thermal_conductivity.model) {
	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		C(0, 0) = C(1, 1) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, time, temp);
		C(0, 1) = C(1, 0) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = _CD(1, 1) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, time, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, time, temp);
		C(0, 1) = C(1, 0) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, time, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, time, temp);
		C(1, 0) = C(0, 1) = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate(p, time, temp);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = derivation(mat->thermal_conductivity.values.get(0, 1), temp / 1e4);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, time, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, time, temp);
		C(0, 1) = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate(p, time, temp);
		C(1, 0) = mat->thermal_conductivity.values.get(1, 0).evaluator->evaluate(p, time, temp);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(0, 1) = derivation(mat->thermal_conductivity.values.get(0, 1), temp / 1e4);
			_CD(1, 0) = derivation(mat->thermal_conductivity.values.get(1, 0), temp / 1e4);
		}
		break;
	default:
		ESINFO(ERROR) << "Advection diffusion 2D not supports set material model";
	}

	if (tangentCorrection) {
		TCDT.multiply(T, _CD * T, 1, 0, true, false);
		CD(node, 0) += phase * TCDT(0, 0);
		CD(node, 1) += phase * TCDT(1, 1);
		CD(node, 2) += phase * TCDT(0, 1);
		CD(node, 3) += phase * TCDT(1, 0);
	}

	TCT.multiply(T, C * T, 1, 0, true, false);
	K(node, 0) += phase * TCT(0, 0);
	K(node, 1) += phase * TCT(1, 1);
	K(node, 2) += phase * TCT(0, 1);
	K(node, 3) += phase * TCT(1, 0);
}

void HeatTransfer2DKernel::processElement(Matrices matrices, const SolverParameters &parameters, const ElementIterator &iterator, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	esint size = iterator.element->nodes;

	const std::vector<DenseMatrix> &N = *(iterator.element->N);
	const std::vector<DenseMatrix> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	bool CAU = info::ecf->heat_transfer_3d.stabilization == HeatTransferConfiguration::STABILIZATION::CAU;
	bool tangentCorrection = parameters.tangentMatrixCorrection;

	DenseMatrix Ce(2, 2), coordinates(size, 2), J(2, 2), invJ(2, 2), dND;
	double detJ, tauK, xi = 1, C1 = 1, C2 = 6;
	DenseMatrix f(size, 1);
	DenseMatrix U(size, 2);
	DenseMatrix m(size, 1);
	DenseMatrix T(size, 1);
	DenseMatrix thickness(size, 1), K(size, 4);
	DenseMatrix gpThickness(1, 1), gpK(1, 4), gpM(1, 1);
	DenseMatrix tangentK, BT, BTN, gpCD, CD, CDBTN, CDe;
	DenseMatrix gKe(size, size);

	if (tangentCorrection) {
		CD.resize(size, 4);
		CDe.resize(2, 2);
	}

	const MaterialBaseConfiguration *phase1, *phase2;
	if (iterator.material->phase_change) {
		phase1 = &iterator.material->phases.find(1)->second;
		phase2 = &iterator.material->phases.find(2)->second;
	}

	for (int n = 0; n < size; n++) {
		T(n, 0) = iterator.temperature[n];
		coordinates(n, 0) = iterator.coordinates[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates[2 * n + 1];
		thickness(n, 0) = iterator.thickness[n];
		if (iterator.material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, iterator.material->phase_change_temperature - iterator.material->transition_interval / 2, iterator.material->phase_change_temperature + iterator.material->transition_interval / 2, T(n, 0), iterator.material->smooth_step_order);
			assembleMaterialMatrix(n, iterator.coordinates + 2 * n, phase1, phase, time::current, T(n, 0), K, CD, tangentCorrection);
			assembleMaterialMatrix(n, iterator.coordinates + 2 * n, phase2, (1 - phase), time::current, T(n, 0), K, CD, tangentCorrection);
			double dens1, dens2, hc1, hc2;
			phase1->density.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &dens1);
			phase2->density.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &dens2);
			phase1->heat_capacity.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &hc1);
			phase2->heat_capacity.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &hc2);

			m(n, 0) = (phase * dens1 + (1 - phase) * dens2) * (phase * hc1 + (1 - phase) * hc2 + iterator.material->latent_heat * derivation) * iterator.thickness[0];
		} else {
			assembleMaterialMatrix(n, iterator.coordinates + 2 * n, iterator.material, 1, time::current, T(n, 0), K, CD, tangentCorrection);
			double dens, hc;
			iterator.material->density.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &dens);
			iterator.material->heat_capacity.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &hc);
			m(n, 0) = dens * hc * thickness(n, 0);
		}

		U(n, 0) = iterator.motion[2 * n + 0] * m(n, 0);
		U(n, 1) = iterator.motion[2 * n + 1] * m(n, 0);
		f(n, 0) = iterator.heat[n];
	}

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if ((matrices & Matrices::K) || ((matrices & Matrices::R) && parameters.timeIntegrationConstantK != 0)) {
		Ke.resize(size, size);
		Ke = 0;
	}
	if ((matrices & Matrices::M) || ((matrices & Matrices::R) && parameters.timeIntegrationConstantM != 0)) {
		Me.resize(size, size);
		Me = 0;
	}
	if (matrices & Matrices::R) {
		Re.resize(size, 1);
		Re = 0;
	}
	if (matrices & Matrices::f) {
		fe.resize(size, 1);
		fe = 0;
	}

	if (tangentCorrection) {
		tangentK.resize(size, size);
	}

	DenseMatrix g(1, 2), u(1, 2), v(1, 2), re(1, size);
	double normGradN = 0;

	if ((matrices & Matrices::M) && info::ecf->heat_transfer_3d.diffusion_split) {
		g(0, 0) = iterator.gradient[0];
		g(0, 1) = iterator.gradient[1];;
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J.values());
		inverse2x2(J.values(), invJ.values(), detJ);

		gpThickness.multiply(N[gp], thickness);
		gpK.multiply(N[gp], K);
		if (tangentCorrection) {
			gpCD.multiply(N[gp], CD);
			CDe(0, 0) = gpCD(0, 0);
			CDe(1, 1) = gpCD(0, 1);
			CDe(0, 1) = gpCD(0, 2);
			CDe(1, 0) = gpCD(0, 3);
		}
		gpM.multiply(N[gp], m);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(0, 1) = gpK(0, 2);
		Ce(1, 0) = gpK(0, 3);

		dND.multiply(invJ, dN[gp]);

		DenseMatrix b_e(1, size), b_e_c(1, size), g_e(1, size);
		b_e.multiply(u, dND, 1, 0);
		g_e.multiply(g, dND, 1, 0);

		if (CAU) {
			normGradN = dND.norm();
			if (normGradN >= 1e-12) {
				for (size_t i = 0; i < re.columns(); i++) {
					re(0, i) = b_e(0, i) - f(i, 0);
				}
				DenseMatrix ReBt(1, 2);
				ReBt.multiply(re, dND, 1 / pow(normGradN, 2), 0, false, true);
				for (size_t i = 0; i < ReBt.columns(); i++) {
					v(0, i) = u(0, i) - ReBt(0, i);
				}
			} else {
				v = u;
			}
		}


		double norm_u_e = u.norm();
		double h_e = 0, tau_e = 0, konst = 0, gh_e = 0;
		double C_e = 0;

		if (info::ecf->heat_transfer_3d.diffusion_split && g.norm() != 0) {
			gh_e = 2 * g.norm() / g_e.norm();
			tauK = (C1 * gh_e * gh_e) / (Ce(0, 0) * C2 + gh_e * gh_e * (gpM(0, 0) / time::shift));
			xi = std::max(1., 1 / (1 - tauK * gpM(0, 0) / time::shift));
		}

		if (norm_u_e != 0) {
			h_e = 2 * norm_u_e / b_e.norm();
			double P_e = h_e * norm_u_e / (2 * Ce(0, 0));
			tau_e = std::max(0.0, 1 - 1 / P_e);
			konst = h_e * tau_e / (2 * norm_u_e);

			if (CAU) {
				DenseMatrix u_v(1, 2);
				u_v(0, 0) = u(0, 0) - v(0, 0);
				u_v(0, 1) = u(0, 1) - v(0, 1);
				b_e_c.multiply(u_v, dND, 1, 0);
				double norm_u_v = u_v.norm();
				double h_e_c = 2 * norm_u_v / b_e_c.norm();
				double P_e_c = h_e_c * norm_u_v / (2 * Ce.norm());
				double tau_e_c = std::max(0.0, 1 - 1 / P_e_c);

				double konst1 = re.norm() / normGradN;
				double konst2 = tau_e * h_e != 0 ? tau_e_c * h_e_c / (tau_e * h_e) : 0;
				if (konst1 / norm_u_e < konst2) {
					C_e = tau_e * h_e * konst1 * (konst2 - konst1 / norm_u_e) / 2;
				} else {
					C_e = 0;
				}
			}
		}

		Ce(0, 0) += info::ecf->heat_transfer_3d.sigma * h_e * norm_u_e;
		Ce(1, 1) += info::ecf->heat_transfer_3d.sigma * h_e * norm_u_e;

		if (matrices & (Matrices::M | Matrices::R)) {
			Me.multiply(N[gp], N[gp], detJ * gpM(0, 0) * weighFactor[gp], 1, true);
		}
		if (matrices & (Matrices::K | Matrices::R)) {
			if (tangentCorrection) {
				BT.multiply(dND, T);
				BTN.multiply(BT, N[gp]);
				CDBTN.multiply(CDe, BTN);
				tangentK.multiply(dND, CDBTN,  detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
			}
			if (info::ecf->heat_transfer_3d.diffusion_split) {
				gKe.multiply(dND, Ce * dND, detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
				gKe.multiply(N[gp], b_e, detJ * weighFactor[gp], 1, true);
				if (konst * weighFactor[gp] * detJ != 0) {
					gKe.multiply(b_e, b_e, konst * weighFactor[gp] * detJ, 1, true);
				}
				if (CAU) {
					gKe.multiply(dND, dND, C_e * weighFactor[gp] * detJ, 1, true);
				}
			}
			Ke.multiply(dND, Ce * dND, xi * detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
			Ke.multiply(N[gp], b_e, xi * detJ * weighFactor[gp], 1, true);
			if (konst * weighFactor[gp] * detJ != 0) {
				Ke.multiply(b_e, b_e, xi * konst * weighFactor[gp] * detJ, 1, true);
			}
			if (CAU) {
				Ke.multiply(dND, dND, xi * C_e * weighFactor[gp] * detJ, 1, true);
			}
		}

		if (matrices & Matrices::f) {
			for (esint i = 0; i < size; i++) {
				fe(i, 0) += detJ * weighFactor[gp] * N[gp](0, i) * f(i, 0);
				if (norm_u_e != 0) {
					fe(i, 0) += detJ * weighFactor[gp] * h_e * tau_e * b_e(0, i) * f(i, 0) / (2 * norm_u_e);
				}
			}
		}
	}

	if ((matrices & Matrices::M) && info::ecf->heat_transfer_3d.diffusion_split) {
		DenseMatrix T1, T2;
		T1.multiply(Ke, T, 1, 0);
		T2.multiply(gKe, T, 1, 0);
		for (esint i = 0; i < size; i++) {
			fe(i, 0) += T1(i, 0) - T2(i, 0);
		}
	}

	if (matrices & Matrices::R) {
		Re.multiply(Ke, T, parameters.timeIntegrationConstantK, 0);
		Re.multiply(Me, T, parameters.timeIntegrationConstantM, 1);
		Re.multiply(Ke, T, 1, 0);
		Re.multiply(Me, T, 1, 1);
		if (!(matrices & Matrices::K)) {
			Ke.resize(0, 0);
		}
		if (!(matrices & Matrices::M)) {
			Me.resize(0, 0);
		}
	}

	if (tangentCorrection) {
		Ke += tangentK;
	}
}

void HeatTransfer2DKernel::processEdge(Matrices matrices, const SolverParameters &parameters, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const
{
	if (!iterator.convection && !iterator.heatFlow && !iterator.heatFlux && !iterator.radiation) {
		Ke.resize(0, 0);
		fe.resize(0, 0);
		return;
	}
	if (!(matrices & (Matrices::K | Matrices::f))) {
		Ke.resize(0, 0);
		fe.resize(0, 0);
		return;
	}

	esint size = iterator.element->nodes;

	const std::vector<DenseMatrix> &N = *(iterator.element->N);
	const std::vector<DenseMatrix> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	DenseMatrix coordinates(size, 2), dND(1, 2), q(size, 1), htc(size, 1), thickness(size, 1), emiss(size, 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpThickness(1, 1), gpEmiss(1, 1);

	Ke.resize(0, 0);
	fe.resize(0, 0);

	if (matrices & Matrices::f) {
		fe.resize(size, 1);
		fe = 0;
	}

	if (iterator.convection || iterator.radiation) {
		Ke.resize(size, size);
		Ke = 0;
	}

	for (int n = 0; n < size; n++) {
		double temp = iterator.temperature[n];
		coordinates(n, 0) = iterator.coordinates[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates[2 * n + 1];
		thickness(n, 0) = iterator.thickness[n];

		if (iterator.convection) {
			double text = iterator.externalTemperature[n];
			htc(n, 0) = iterator.htc[n];

			if (time::iteration) {
				q(n, 0) += htc(n, 0) * (text - temp);
			} else {
				q(n, 0) += htc(n, 0) * (text);
			}
		}

		if (iterator.radiation) {
			emiss(n, 0) = CONST_Stefan_Boltzmann * iterator.emissivity[n];
			q(n, 0) += emiss(n, 0) * (pow(iterator.externalTemperature[n], 4) - pow(temp, 4));
			emiss(n, 0) *= 4 * temp * temp * temp;
		}
		if (iterator.heatFlow) {
			q(n, 0) += iterator.heatFlow[n] / iterator.regionArea;
		}
		if (iterator.heatFlux) {
			q(n, 0) += iterator.heatFlux[n];
		}

		q(n, 0) *= thickness(n, 0);
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);

		if (iterator.convection || iterator.radiation) {
			gpThickness.multiply(N[gp], thickness);
		}
		if (iterator.convection) {
			gpHtc.multiply(N[gp], htc);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
		}
		if (iterator.radiation) {
			gpEmiss.multiply(N[gp], emiss);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0), 1, true);
		}
		for (esint i = 0; i < size; i++) {
			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % size) * gpQ(0, 0);
		}
	}
}

void HeatTransfer2DKernel::processSolution(const SolutionIterator &iterator)
{
	esint size = iterator.element->nodes;

	const std::vector<DenseMatrix> &N = *(iterator.element->N);
	const std::vector<DenseMatrix> &dN = *(iterator.element->dN);

	DenseMatrix Ce(2, 2), coordinates(size, 2), J(2, 2), invJ(2, 2), dND, T(size, 1);
	double detJ, m, norm_u_e, h_e;
	DenseMatrix thickness(size, 1), U(size, 2), K(size, 4), gpK(1, 4), CD;
	DenseMatrix u(1, 2), matFlux(2, 1), matGradient(2, 1);

	const MaterialBaseConfiguration *phase1, *phase2;
	if (iterator.material->phase_change) {
		phase1 = &iterator.material->phases.find(1)->second;
		phase2 = &iterator.material->phases.find(2)->second;
	}

	for (int n = 0; n < size; n++) {
		T(n, 0) = iterator.temperature[n];
		coordinates(n, 0) = iterator.coordinates[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates[2 * n + 1];
		thickness(n, 0) = iterator.thickness[n];
		if (iterator.material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, iterator.material->phase_change_temperature - iterator.material->transition_interval / 2, iterator.material->phase_change_temperature + iterator.material->transition_interval / 2, T(n, 0), iterator.material->smooth_step_order);
			assembleMaterialMatrix(n, iterator.coordinates + 2 * n, phase1, phase, time::current, T(n, 0), K, CD, false);
			assembleMaterialMatrix(n, iterator.coordinates + 2 * n, phase2, (1 - phase), time::current, T(n, 0), K, CD, false);
			double dens1, dens2, hc1, hc2;
			phase1->density.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &dens1);
			phase2->density.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &dens2);
			phase1->heat_capacity.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &hc1);
			phase2->heat_capacity.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &hc2);

			m = (phase * dens1 + (1 - phase) * dens2) * (phase * hc1 + (1 - phase) * hc2 + iterator.material->latent_heat * derivation) * iterator.thickness[0];
			if (iterator.material->phase_change) {
				*iterator.phase = phase;
				*iterator.latentHeat = iterator.material->latent_heat * derivation;
			}
		} else {
			assembleMaterialMatrix(n, iterator.coordinates + 2 * n, iterator.material, 1, time::current, T(n, 0), K, CD, false);
			double dens, hc;
			iterator.material->density.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &dens);
			iterator.material->heat_capacity.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &hc);
			m = dens * hc * thickness(n, 0);
		}

		U(n, 0) = iterator.motion[2 * n + 0] * m;
		U(n, 1) = iterator.motion[2 * n + 1] * m;
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J.values());
		inverse2x2(J.values(), invJ.values(), detJ);

		gpK.multiply(N[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(0, 1) = gpK(0, 2);
		Ce(1, 0) = gpK(0, 3);

		dND.multiply(invJ, dN[gp]);

		norm_u_e = u.norm();
		h_e = 0;

		if (norm_u_e != 0) {
			DenseMatrix b_e(1, size);
			b_e.multiply(u, dND, 1, 0);
			h_e = 2 * norm_u_e / b_e.norm();
		}

		Ce(0, 0) += info::ecf->heat_transfer_3d.sigma * h_e * norm_u_e;
		Ce(1, 1) += info::ecf->heat_transfer_3d.sigma * h_e * norm_u_e;

		if (info::ecf->output.results_selection.gradient) {
			matGradient.multiply(dND, T, 1, 1);
		}
		if (info::ecf->output.results_selection.flux) {
			matFlux.multiply(Ce, dND * T, 1, 1);
		}
	}

	if (info::ecf->output.results_selection.gradient) {
		*(iterator.gradient + 0) = matGradient(0, 0) / N.size();
		*(iterator.gradient + 1) = matGradient(1, 0) / N.size();
	}

	if (info::ecf->output.results_selection.flux) {
		*(iterator.flux + 0) = matFlux(0, 0) / N.size();
		*(iterator.flux+ 1) = matFlux(1, 0) / N.size();
	}
}





