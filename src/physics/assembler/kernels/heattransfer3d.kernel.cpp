
#include "heattransfer3d.kernel.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/matrices/denseMatrix.h"
#include "../../../basis/evaluator/evaluator.h"

#include "../../../config/ecf/physics/heattransfer.h"
#include "../../../globals/time.h"

#include "../../../mesh/elements/element.h"
#include "../../dataholder.h"

using namespace espreso;

// TODO: create file with constants
#define CONST_Stefan_Boltzmann 5.6703e-8

using namespace espreso;

HeatTransfer3DKernel::HeatTransfer3DKernel(const HeatTransferGlobalSettings &settings, const HeatTransferOutputSettings &output)
: HeatTransferKernel(settings, output)
{

}

void HeatTransfer3DKernel::assembleMaterialMatrix(eslocal node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const
{
	auto d2r = [] (double degree) -> double {
		return M_PI * degree / 180;
	};

	Point p(coordinates[0], coordinates[1], coordinates[2]);

	Point sin, cos;
	switch (mat->coordinate_system.type) {
	case CoordinateSystemConfiguration::TYPE::CARTESIAN:

		cos.x = std::cos(d2r(mat->coordinate_system.rotation.x.evaluator->evaluate(p, time, temp)));
		cos.y = std::cos(d2r(mat->coordinate_system.rotation.y.evaluator->evaluate(p, time, temp)));
		cos.z = std::cos(d2r(mat->coordinate_system.rotation.z.evaluator->evaluate(p, time, temp)));

		sin.x = std::sin(d2r(mat->coordinate_system.rotation.x.evaluator->evaluate(p, time, temp)));
		sin.y = std::sin(d2r(mat->coordinate_system.rotation.y.evaluator->evaluate(p, time, temp)));
		sin.z = std::sin(d2r(mat->coordinate_system.rotation.z.evaluator->evaluate(p, time, temp)));

		break;

	case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: {

		Point origin(
				mat->coordinate_system.center.x.evaluator->evaluate(p, time, temp),
				mat->coordinate_system.center.y.evaluator->evaluate(p, time, temp),
				0);

		double rotation = std::atan2((p.y - origin.y), (p.x - origin.x));

		cos.x = 1.0;
		cos.y = 1.0;
		cos.z = std::cos(rotation);

		sin.x = 0.0;
		sin.y = 0.0;
		sin.z = std::sin(rotation);

	} break;

	case CoordinateSystemConfiguration::TYPE::SPHERICAL: {

		Point origin(
				mat->coordinate_system.center.x.evaluator->evaluate(p, time, temp),
				mat->coordinate_system.center.y.evaluator->evaluate(p, time, temp),
				mat->coordinate_system.center.z.evaluator->evaluate(p, time, temp));

		double azimut = std::atan2((p.y - origin.y), (p.x - origin.x));
		double r = std::sqrt(pow((p.x - origin.x), 2) + pow((p.y - origin.y), 2) + pow((p.z - origin.z), 2));
		double elevation = 0.0;

		if (r < 1e-12) {
			elevation = 0.0;
		} else {
			elevation = std::atan2(std::sqrt(pow((p.z - origin.z), 2) + pow((p.x - origin.x), 2)), (p.y - origin.y));
		}

		cos.x = 1.0;
		cos.y = std::cos(elevation);
		cos.z = std::cos(azimut);

		sin.x = 0.0;
		sin.y = std::sin(elevation);
		sin.z = std::sin(azimut);

	} break;

	}

	DenseMatrix TCT(3, 3), T(3, 3), C(3, 3), _CD, TCDT;

	T(0, 0) = cos.y * cos.z;                         T(0, 1) = cos.y * sin.z;                         T(0, 2) = -sin.y;
	T(1, 0) = cos.z * sin.x * sin.y - cos.x * sin.z; T(1, 1) = cos.x * cos.z + sin.x * sin.y * sin.z; T(1, 2) = cos.y * sin.x;
	T(2, 0) = sin.x * sin.z + cos.x * cos.z * sin.y; T(2, 1) = cos.x * sin.y * sin.z - cos.z * sin.x; T(2, 2) = cos.x * cos.y;

	if (tangentCorrection) {
		_CD.resize(3, 3);
		TCDT.resize(3, 3);
	}

	auto derivation = [&] (const ECFExpression &expression, double h) {
		return (
				expression.evaluator->evaluate(p, time, temp + h) -
				expression.evaluator->evaluate(p, time, temp - h)
				) / (2 * h);
	};

	switch (mat->thermal_conductivity.model) {
	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		C(0, 0) = C(1, 1) = C(2, 2) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, time, temp);
		C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = _CD(1, 1) = _CD(2, 2) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(0, 1) = _CD(0, 2) = _CD(1, 0) = _CD(1, 2) = _CD(2, 0) = _CD(2, 1) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, time, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, time, temp);
		C(2, 2) = mat->thermal_conductivity.values.get(2, 2).evaluator->evaluate(p, time, temp);
		C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(2, 2) = derivation(mat->thermal_conductivity.values.get(2, 2), temp / 1e4);
			_CD(0, 1) = _CD(0, 2) = _CD(1, 0) = _CD(1, 2) = _CD(2, 0) = _CD(2, 1) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, time, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, time, temp);
		C(2, 2) = mat->thermal_conductivity.values.get(2, 2).evaluator->evaluate(p, time, temp);
		C(0, 1) = C(1, 0) = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate(p, time, temp);
		C(0, 2) = C(2, 0) = mat->thermal_conductivity.values.get(0, 2).evaluator->evaluate(p, time, temp);
		C(1, 2) = C(2, 1) = mat->thermal_conductivity.values.get(1, 2).evaluator->evaluate(p, time, temp);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(2, 2) = derivation(mat->thermal_conductivity.values.get(2, 2), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = derivation(mat->thermal_conductivity.values.get(0, 1), temp / 1e4);
			_CD(0, 2) = _CD(2, 0) = derivation(mat->thermal_conductivity.values.get(0, 2), temp / 1e4);
			_CD(1, 2) = _CD(2, 1) = derivation(mat->thermal_conductivity.values.get(1, 2), temp / 1e4);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, time, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, time, temp);
		C(2, 2) = mat->thermal_conductivity.values.get(2, 2).evaluator->evaluate(p, time, temp);
		C(0, 1) = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate(p, time, temp);
		C(0, 2) = mat->thermal_conductivity.values.get(0, 2).evaluator->evaluate(p, time, temp);
		C(1, 0) = mat->thermal_conductivity.values.get(1, 0).evaluator->evaluate(p, time, temp);
		C(1, 2) = mat->thermal_conductivity.values.get(1, 2).evaluator->evaluate(p, time, temp);
		C(2, 0) = mat->thermal_conductivity.values.get(2, 0).evaluator->evaluate(p, time, temp);
		C(2, 1) = mat->thermal_conductivity.values.get(2, 1).evaluator->evaluate(p, time, temp);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(2, 2) = derivation(mat->thermal_conductivity.values.get(2, 2), temp / 1e4);
			_CD(0, 1) = derivation(mat->thermal_conductivity.values.get(0, 1), temp / 1e4);
			_CD(0, 2) = derivation(mat->thermal_conductivity.values.get(0, 2), temp / 1e4);
			_CD(1, 2) = derivation(mat->thermal_conductivity.values.get(1, 2), temp / 1e4);
			_CD(1, 0) = derivation(mat->thermal_conductivity.values.get(1, 0), temp / 1e4);
			_CD(2, 0) = derivation(mat->thermal_conductivity.values.get(2, 0), temp / 1e4);
			_CD(2, 1) = derivation(mat->thermal_conductivity.values.get(2, 1), temp / 1e4);
		}
		break;
	default:
		ESINFO(ERROR) << "Advection diffusion 3D not supports set material model";
	}

	TCT.multiply(T, C * T, 1, 0, true, false);
	if (tangentCorrection) {
		TCDT.multiply(T, _CD * T, 1, 0, true, false);
		CD(node, 0) += phase * TCDT(0, 0);
		CD(node, 1) += phase * TCDT(1, 1);
		CD(node, 2) += phase * TCDT(2, 2);
		CD(node, 3) += phase * TCDT(0, 1);
		CD(node, 4) += phase * TCDT(0, 2);
		CD(node, 5) += phase * TCDT(1, 0);
		CD(node, 6) += phase * TCDT(1, 2);
		CD(node, 7) += phase * TCDT(2, 0);
		CD(node, 8) += phase * TCDT(2, 1);
	}

	K(node, 0) += phase * TCT(0, 0);
	K(node, 1) += phase * TCT(1, 1);
	K(node, 2) += phase * TCT(2, 2);
	K(node, 3) += phase * TCT(0, 1);
	K(node, 4) += phase * TCT(0, 2);
	K(node, 5) += phase * TCT(1, 0);
	K(node, 6) += phase * TCT(1, 2);
	K(node, 7) += phase * TCT(2, 0);
	K(node, 8) += phase * TCT(2, 1);
}

void HeatTransfer3DKernel::processElement(Matrices matrices, const ElementIterator &iterator, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	eslocal size = iterator.element->nodes;

	const std::vector<DenseMatrix> &N = *(iterator.element->N);
	const std::vector<DenseMatrix> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	bool CAU = _settings.stabilization == HeatTransferConfiguration::STABILIZATION::CAU;
	bool tangentCorrection = (matrices & Matrices::K); // && step.tangentMatrixCorrection;

	DenseMatrix Ce(3, 3), coordinates(size, 3), J(3, 3), invJ(3, 3), dND;
	double detJ, temp, tauK, xi = 1, C1 = 1, C2 = 6;
	DenseMatrix f(size, 1);
	DenseMatrix U(size, 3);
	DenseMatrix m(size, 1);
	DenseMatrix T(size, 1);
	DenseMatrix K(size, 9);
	DenseMatrix gpK(1, 9), gpM(1, 1);
	DenseMatrix tangentK, BT, BTN, gpCD, CD, CDBTN, CDe;
	DenseMatrix gKe(size, size);

	if (tangentCorrection) {
		CD.resize(size, 9);
		CDe.resize(3, 3);
	}

	const MaterialBaseConfiguration *phase1, *phase2;
	if (iterator.material->phase_change) {
		phase1 = &iterator.material->phases.find(1)->second;
		phase2 = &iterator.material->phases.find(2)->second;
	}

	for (int n = 0; n < size; n++) {
		T(n, 0) = iterator.temperature[n];
		coordinates(n, 0) = iterator.coordinates[3 * n + 0];
		coordinates(n, 1) = iterator.coordinates[3 * n + 1];
		coordinates(n, 2) = iterator.coordinates[3 * n + 2];
		if (iterator.material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, iterator.material->phase_change_temperature - iterator.material->transition_interval / 2, iterator.material->phase_change_temperature + iterator.material->transition_interval / 2, T(n, 0), iterator.material->smooth_step_order);
			assembleMaterialMatrix(n, iterator.coordinates + 3 * n, phase1, phase, time::current, T(n, 0), K, CD, tangentCorrection);
			assembleMaterialMatrix(n, iterator.coordinates + 3 * n, phase2, (1 - phase), time::current, T(n, 0), K, CD, tangentCorrection);
			double dens1, dens2, hc1, hc2;
			phase1->density.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &dens1);
			phase2->density.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &dens2);
			phase1->heat_capacity.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &hc1);
			phase2->heat_capacity.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &hc2);

			m(n, 0) = (phase * dens1 + (1 - phase) * dens2) * (phase * hc1 + (1 - phase) * hc2 + iterator.material->latent_heat * derivation);
		} else {
			assembleMaterialMatrix(n, iterator.coordinates + 3 * n, iterator.material, 1, time::current, T(n, 0), K, CD, tangentCorrection);
			double dens, hc;
			iterator.material->density.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &dens);
			iterator.material->heat_capacity.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &hc);
			m(n, 0) = dens * hc;
		}

		U(n, 0) = iterator.motion[3 * n + 0] * m(n, 0);
		U(n, 1) = iterator.motion[3 * n + 1] * m(n, 0);
		U(n, 2) = iterator.motion[3 * n + 2] * m(n, 0);
		f(n, 0) = iterator.heat[n];
	}

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if ((matrices & Matrices::K) || ((matrices & Matrices::R))) { // && step.timeIntegrationConstantK != 0)) {
		Ke.resize(size, size);
		Ke = 0;
	}
	if ((matrices & Matrices::M) || ((matrices & Matrices::R))) { // && step.timeIntegrationConstantM != 0)) {
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

	DenseMatrix g(1, 3), u(1, 3), v(1, 3), re(1, size);
	double normGradN = 0;

	if ((matrices & Matrices::M) && _settings.diffusion_split) {
		g(0, 0) = iterator.gradient[0];
		g(0, 1) = iterator.gradient[1];
		g(0, 2) = iterator.gradient[2];
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J.values());
		inverse3x3(J.values(), invJ.values(), detJ);

		gpK.multiply(N[gp], K);
		if (tangentCorrection) {
			gpCD.multiply(N[gp], CD);
			CDe(0, 0) = gpCD(0, 0);
			CDe(1, 1) = gpCD(0, 1);
			CDe(2, 2) = gpCD(0, 2);
			CDe(0, 1) = gpCD(0, 3);
			CDe(0, 2) = gpCD(0, 4);
			CDe(1, 0) = gpCD(0, 5);
			CDe(1, 2) = gpCD(0, 6);
			CDe(2, 0) = gpCD(0, 7);
			CDe(2, 1) = gpCD(0, 8);
		}
		gpM.multiply(N[gp], m);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(2, 2) = gpK(0, 2);
		Ce(0, 1) = gpK(0, 3);
		Ce(0, 2) = gpK(0, 4);
		Ce(1, 2) = gpK(0, 6);
		Ce(1, 0) = gpK(0, 5);
		Ce(2, 0) = gpK(0, 7);
		Ce(2, 1) = gpK(0, 8);

		dND.multiply(invJ, dN[gp]);

		DenseMatrix b_e(1, size), b_e_c(1, size), g_e(1, size);
		b_e.multiply(u, dND, 1, 0);
		g_e.multiply(g, dND, 1, 0);

		if (CAU) {
			normGradN = dND.norm();
			if (normGradN >= 1e-12) {
				for (size_t i = 0; i < re.columns(); i++) {
					re(0, i) = b_e(0, i) - f(0, i);
				}
				DenseMatrix ReBt(1, 3);
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

		if ((matrices & Matrices::M) && _settings.diffusion_split && g.norm() != 0) {
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
				DenseMatrix u_v(1, 3);
				u_v(0, 0) = u(0, 0) - v(0, 0);
				u_v(0, 1) = u(0, 1) - v(0, 1);
				u_v(0, 2) = u(0, 2) - v(0, 2);
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

		Ce(0, 0) += _settings.sigma * h_e * norm_u_e;
		Ce(1, 1) += _settings.sigma * h_e * norm_u_e;
		Ce(2, 2) += _settings.sigma * h_e * norm_u_e;

		if (matrices & (Matrices::M | Matrices::R)) {
			Me.multiply(N[gp], N[gp], detJ * gpM(0, 0) * weighFactor[gp], 1, true);
		}
		if (matrices & (Matrices::K | Matrices::R)) {
			if (tangentCorrection) {
				BT.multiply(dND, T);
				BTN.multiply(BT, N[gp]);
				CDBTN.multiply(CDe, BTN);
				tangentK.multiply(dND, CDBTN,  detJ * weighFactor[gp], 1, true);
			}
			if ((matrices & Matrices::M) && _settings.diffusion_split) {
				gKe.multiply(dND, Ce * dND, detJ * weighFactor[gp], 1, true);
				gKe.multiply(N[gp], b_e, detJ * weighFactor[gp], 1, true);
				if (konst * weighFactor[gp] * detJ != 0) {
					gKe.multiply(b_e, b_e, konst * weighFactor[gp] * detJ, 1, true);
				}
				if (CAU) {
					gKe.multiply(dND, dND, C_e * weighFactor[gp] * detJ, 1, true);
				}
			}
			Ke.multiply(dND, Ce * dND, xi * detJ * weighFactor[gp], 1, true);
			Ke.multiply(N[gp], b_e, xi * detJ * weighFactor[gp], 1, true);
			if (konst * weighFactor[gp] * detJ != 0) {
				Ke.multiply(b_e, b_e, xi * konst * weighFactor[gp] * detJ, 1, true);
			}
			if (CAU) {
				Ke.multiply(dND, dND, xi * C_e * weighFactor[gp] * detJ, 1, true);
			}
		}

		if (matrices & Matrices::f) {
			for (eslocal i = 0; i < size; i++) {
				fe(i, 0) += detJ * weighFactor[gp] * N[gp](0, i) * f(i, 0);
				if (norm_u_e != 0) {
					fe(i, 0) += detJ * weighFactor[gp] * h_e * tau_e * b_e(0, i) * f(i, 0) / (2 * norm_u_e);
				}
			}
		}
	}

	if ((matrices & Matrices::M) && _settings.diffusion_split) {
		DenseMatrix T1, T2;
		T1.multiply(Ke, T, 1, 0);
		T2.multiply(gKe, T, 1, 0);
		for (eslocal i = 0; i < size; i++) {
			fe(i, 0) += T1(i, 0) - T2(i, 0);
		}
	}

	if (matrices & Matrices::R) {
//		Re.multiply(Ke, T, step.timeIntegrationConstantK, 0);
//		Re.multiply(Me, T, step.timeIntegrationConstantM, 1);
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

void HeatTransfer3DKernel::processFace(Matrices matrices, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const
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

	eslocal size = iterator.element->nodes;

	const std::vector<DenseMatrix> &N = *(iterator.element->N);
	const std::vector<DenseMatrix> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	DenseMatrix coordinates(size, 3), dND(1, 3), q(size, 1), htc(size, 1), flow(size, 1), emiss(size, 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpFlow(1, 1), gpEmiss(1, 1);

	double area = iterator.regionArea, temp, text;
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

	for (size_t n = 0; n < size; n++) {
		double temp = iterator.temperature[n];
		coordinates(n, 0) = iterator.coordinates[3 * n + 0];
		coordinates(n, 1) = iterator.coordinates[3 * n + 1];
		coordinates(n, 2) = iterator.coordinates[3 * n + 2];
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
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		double J = va.norm();

		gpQ.multiply(N[gp], q);
		if (iterator.convection) {
			gpHtc.multiply(N[gp], htc);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
		}
		if (iterator.radiation) {
			gpEmiss.multiply(N[gp], emiss);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0), 1, true);
		}
		for (eslocal i = 0; i < size; i++) {
			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % size) * gpQ(0, 0);
		}
	}
}

void HeatTransfer3DKernel::processEdge(Matrices matrices, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const
{
//	if (!(e->hasProperty(Property::EXTERNAL_TEMPERATURE, _step->step) ||
//		e->hasProperty(Property::HEAT_FLOW, _step->step) ||
//		e->hasProperty(Property::HEAT_FLUX, _step->step))) {
//
//		Ke.resize(0, 0);
//		Me.resize(0, 0);
//		Re.resize(0, 0);
//		fe.resize(0, 0);
//		return;
//	}
//	if (!(matrices & (Matrices::K | Matrices::f))) {
//		Ke.resize(0, 0);
//		Me.resize(0, 0);
//		Re.resize(0, 0);
//		fe.resize(0, 0);
//		return;
//	}
//
//	DenseMatrix coordinates(e->nodes(), 3), dND(1, 3), q(e->nodes(), 1), htc(e->nodes(), 1), flow(e->nodes(), 1), emiss(e->nodes(), 1);
//	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpFlow(1, 1), gpEmiss(1, 1);
//
//	double area = 1, temp;
//	eslocal Ksize = e->nodes();
//	Ke.resize(0, 0);
//	Me.resize(0, 0);
//	Re.resize(0, 0);
//	fe.resize(0, 0);
//
//	if (matrices & Matrices::f) {
//		fe.resize(Ksize, 1);
//		fe = 0;
//	}
//
//	for (size_t r = 0; r < e->regions().size(); r++) {
//		if (_step->step < e->regions()[r]->settings.size() && e->regions()[r]->settings[_step->step].count(Property::HEAT_FLOW)) {
//			area = e->regions()[r]->area;
//			break;
//		}
//	}
//	if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, _step->step)) {
//		Ke.resize(Ksize, Ksize);
//		Ke = 0;
//	}
//
//	const std::vector<DenseMatrix> &dN = e->dN();
//	const std::vector<DenseMatrix> &N = e->N();
//	const std::vector<double> &weighFactor = e->weighFactor();
//
//	const ConvectionConfiguration *convection = NULL;
//	for (size_t r = 0; convection == NULL && r < e->regions().size(); r++) {
//		auto regionit = _configuration.load_steps_settings.at(_step->step + 1).convection.find(e->regions()[r]->name);
//		if (regionit != _configuration.load_steps_settings.at(_step->step + 1).convection.end()) {
//			convection = &regionit->second;
//		}
//	}

//	for (size_t n = 0; n < e->nodes(); n++) {
//		coordinates(n, 0) = _mesh->coordinates()[e->node(n)].x;
//		coordinates(n, 1) = _mesh->coordinates()[e->node(n)].y;
//		coordinates(n, 2) = _mesh->coordinates()[e->node(n)].z;
//
//		temp = solution[offset + SolutionIndex::TEMPERATURE]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(n), e->domains().front()));
//		htc(n, 0) = convection != NULL ? computeHTC(*convection, e, _mesh->coordinates()[e->node(n)], step, temp) : 0;
//
//		if (_step->iteration) {
//			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, _step->step, _mesh->coordinates()[e->node(n)], _step->currentTime, temp, 0) - temp);
//		} else {
//			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, _step->step, _mesh->coordinates()[e->node(n)], _step->currentTime, temp, 0));
//		}
//
//		emiss(n, 0) = CONST_Stefan_Boltzmann * e->getProperty(Property::EMISSIVITY, _step->step, _mesh->coordinates()[e->node(n)], _step->currentTime, temp, 0);
//		q(n, 0) += emiss(n, 0) * (pow(e->getProperty(Property::EXTERNAL_TEMPERATURE, _step->step, _mesh->coordinates()[e->node(n)], _step->currentTime, temp, 0), 4) - pow(temp, 4));
//		q(n, 0) += e->getProperty(Property::HEAT_FLOW, _step->step, _mesh->coordinates()[e->node(n)], _step->currentTime, temp, 0) / area;
//		q(n, 0) += e->getProperty(Property::HEAT_FLUX, _step->step, _mesh->coordinates()[e->node(n)], _step->currentTime, temp, 0);
//
//		emiss(n, 0) *= 4 * temp * temp * temp;
//	}
//
//	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
//		dND.multiply(dN[gp], coordinates);
//		double J = dND.norm();
//		gpQ.multiply(N[gp], q);
//		if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, _step->step)) {
//			gpHtc.multiply(N[gp], htc);
//			gpEmiss.multiply(N[gp], emiss);
//
//			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
//			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0), 1, true);
//		}
//		for (eslocal i = 0; i < Ksize; i++) {
//			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % e->nodes()) * gpQ(0, 0);
//		}
//	}
}

void HeatTransfer3DKernel::processSolution(const SolutionIterator &iterator)
{
	eslocal size = iterator.element->nodes;

	const std::vector<DenseMatrix> &N = *(iterator.element->N);
	const std::vector<DenseMatrix> &dN = *(iterator.element->dN);

	DenseMatrix Ce(3, 3), coordinates(size, 3), J(3, 3), invJ(3, 3), dND, T(size, 1);
	double detJ, norm_u_e, h_e;
	DenseMatrix U(size, 3), K(size, 9), gpK(1, 9), CD;
	DenseMatrix u(1, 3), matFlux(3, 1), matGradient(3, 1);

	const MaterialBaseConfiguration *phase1, *phase2;
	if (iterator.material->phase_change) {
		phase1 = &iterator.material->phases.find(1)->second;
		phase2 = &iterator.material->phases.find(2)->second;
	}

	for (int n = 0; n < size; n++) {
		T(n, 0) = iterator.temperature[n];
		coordinates(n, 0) = iterator.coordinates[3 * n + 0];
		coordinates(n, 1) = iterator.coordinates[3 * n + 1];
		coordinates(n, 2) = iterator.coordinates[3 * n + 2];
		if (iterator.material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, iterator.material->phase_change_temperature - iterator.material->transition_interval / 2, iterator.material->phase_change_temperature + iterator.material->transition_interval / 2, T(n, 0), iterator.material->smooth_step_order);
			assembleMaterialMatrix(n, iterator.coordinates + 3 * n, phase1, phase, time::current, T(n, 0), K, CD, false);
			assembleMaterialMatrix(n, iterator.coordinates + 3 * n, phase2, (1 - phase), time::current, T(n, 0), K, CD, false);
			double dens1, dens2, hc1, hc2;
			phase1->density.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &dens1);
			phase2->density.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &dens2);
			phase1->heat_capacity.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &hc1);
			phase2->heat_capacity.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &hc2);
			if (iterator.material->phase_change) {
				*iterator.phase = phase;
				*iterator.latentHeat = iterator.material->latent_heat * derivation;
			}
		} else {
			assembleMaterialMatrix(n, iterator.coordinates + 3 * n, iterator.material, 1, time::current, T(n, 0), K, CD, false);
		}

		U(n, 0) = iterator.motion[3 * n + 0];
		U(n, 1) = iterator.motion[3 * n + 1];
		U(n, 2) = iterator.motion[3 * n + 2];
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J.values());
		inverse3x3(J.values(), invJ.values(), detJ);

		gpK.multiply(N[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(2, 2) = gpK(0, 2);
		Ce(0, 1) = gpK(0, 3);
		Ce(0, 2) = gpK(0, 4);
		Ce(1, 2) = gpK(0, 6);
		Ce(1, 0) = gpK(0, 5);
		Ce(2, 0) = gpK(0, 7);
		Ce(2, 1) = gpK(0, 8);

		dND.multiply(invJ, dN[gp]);

		norm_u_e = u.norm();
		h_e = 0;

		if (norm_u_e != 0) {
			DenseMatrix b_e(1, size);
			b_e.multiply(u, dND, 1, 0);
			h_e = 2 * norm_u_e / b_e.norm();
		}

		Ce(0, 0) += _settings.sigma * h_e * norm_u_e;
		Ce(1, 1) += _settings.sigma * h_e * norm_u_e;
		Ce(2, 2) += _settings.sigma * h_e * norm_u_e;

		if (_output.gradient) {
			matGradient.multiply(dND, T, 1, 1);
		}
		if (_output.flux) {
			matFlux.multiply(Ce, dND * T, 1, 1);
		}
	}

	if (_output.gradient) {
		*(iterator.gradient + 0) = matGradient(0, 0) / N.size();
		*(iterator.gradient + 1) = matGradient(1, 0) / N.size();
		*(iterator.gradient + 2) = matGradient(2, 0) / N.size();
	}

	if (_output.flux) {
		*(iterator.flux + 0) = matFlux(0, 0) / N.size();
		*(iterator.flux + 1) = matFlux(1, 0) / N.size();
		*(iterator.flux + 2) = matFlux(2, 0) / N.size();
	}
}

