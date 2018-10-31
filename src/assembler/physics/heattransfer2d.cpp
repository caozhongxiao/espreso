
#include "heattransfer2d.h"

#include "../step.h"
#include "../instance.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/evaluator/evaluator.h"
#include "../../basis/matrices/denseMatrix.h"

#include "../../config/ecf/physics/heattransfer.h"
#include "../../config/ecf/output.h"

#include "../../mesh/mesh.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/nodestore.h"
#include "../../mesh/store/boundaryregionstore.h"
#include "../../mesh/store/elementsregionstore.h"

#include "../../solver/generic/SparseMatrix.h"


// TODO: create file with constants
#define CONST_Stefan_Boltzmann 5.6703e-8

using namespace espreso;

HeatTransfer2D::HeatTransfer2D(Mesh *mesh, Instance *instance, Step *step, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration)
: Physics("HEAT TRANSFER 2D", mesh, instance, step, &configuration, 1), HeatTransfer(configuration, propertiesConfiguration)
{
	if (_propertiesConfiguration.gradient || _configuration.diffusion_split) {
		_gradient = _mesh->elements->appendData(2, { "GRADIENT", "GRADIENT_X", "GRADIENT_Y" });
	}

	if (_propertiesConfiguration.flux) {
		_flux = _mesh->elements->appendData(2, { "FLUX", "FLUX_X", "FLUX_Y" });
	}

	bool translationMotion = false;
	for (auto ls = _configuration.load_steps_settings.begin(); ls != _configuration.load_steps_settings.end(); ++ls) {
		if (ls->second.translation_motions.size()) {
			translationMotion = true;
		}
	}
	if (translationMotion) {
		_translationMotion = _mesh->elements->appendData(2, { "TRANSLATION_MOTION", "TRANSLATION_MOTION_X", "TRANSLATION_MOTION_Y" });
	}
}

void HeatTransfer2D::processBEM(eslocal domain, Matrices matrices)
{
	ESINFO(ERROR) << "BEM for 2D heat transfer is not implemented.";
}

void HeatTransfer2D::assembleMaterialMatrix(eslocal node, const Point &p, const MaterialBaseConfiguration *mat, double phase, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const
{
	auto d2r = [] (double degree) -> double {
		return M_PI * degree / 180;
	};

	double cos, sin;
	switch (mat->coordinate_system.type) {
	case CoordinateSystemConfiguration::TYPE::CARTESIAN:
		cos = std::cos(d2r(mat->coordinate_system.rotation.z.evaluator->evaluate(p, _step->currentTime, temp)));
		sin = std::sin(d2r(mat->coordinate_system.rotation.z.evaluator->evaluate(p, _step->currentTime, temp)));
		break;
	case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: {
		Point origin(
				mat->coordinate_system.center.x.evaluator->evaluate(p, _step->currentTime, temp),
				mat->coordinate_system.center.y.evaluator->evaluate(p, _step->currentTime, temp),
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
				expression.evaluator->evaluate(p, _step->currentTime, temp + h) -
				expression.evaluator->evaluate(p, _step->currentTime, temp - h)
				) / (2 * h);
	};

	switch (mat->thermal_conductivity.model) {
	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		C(0, 0) = C(1, 1) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(0, 1) = C(1, 0) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = _CD(1, 1) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		C(0, 1) = C(1, 0) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 0) = C(0, 1) = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate(p, _step->currentTime, temp);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = derivation(mat->thermal_conductivity.values.get(0, 1), temp / 1e4);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		C(0, 1) = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 0) = mat->thermal_conductivity.values.get(1, 0).evaluator->evaluate(p, _step->currentTime, temp);
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

void HeatTransfer2D::processElement(eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	auto nodes = _mesh->elements->procNodes->cbegin() + eindex;
	auto epointer = _mesh->elements->epointers->datatarray()[eindex];
	const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[domain];
	const ECFExpressionVector *translation_motion = NULL;
	Evaluator *heat_source = NULL;
	Evaluator *thick = NULL;
	for (auto it = _configuration.load_steps_settings.at(_step->step + 1).translation_motions.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).translation_motions.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			translation_motion = &it->second;
			break;
		}
	}
	for (auto it = _configuration.load_steps_settings.at(_step->step + 1).heat_source.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).heat_source.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			heat_source = it->second.evaluator;
			break;
		}
	}
	for (auto it = _configuration.thickness.begin(); it != _configuration.thickness.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			thick = it->second.evaluator;
			break;
		}
	}

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	bool CAU = _configuration.stabilization == HeatTransferConfiguration::STABILIZATION::CAU;
	bool tangentCorrection = (matrices & Matrices::K) && _step->tangentMatrixCorrection;

	DenseMatrix Ce(2, 2), coordinates(nodes->size(), 2), J(2, 2), invJ(2, 2), dND;
	double detJ, tauK, xi = 1, C1 = 1, C2 = 6;
	Point center;
	DenseMatrix f(nodes->size(), 1);
	DenseMatrix U(nodes->size(), 2);
	DenseMatrix m(nodes->size(), 1);
	DenseMatrix T(nodes->size(), 1);
	DenseMatrix thickness(nodes->size(), 1), K(nodes->size(), 4);
	DenseMatrix gpThickness(1, 1), gpK(1, 4), gpM(1, 1);
	DenseMatrix tangentK, BT, BTN, gpCD, CD, CDBTN, CDe;
	DenseMatrix gKe(nodes->size(), nodes->size());

	const MaterialConfiguration* material = _mesh->materials[_mesh->elements->material->datatarray()[eindex]];

	if (tangentCorrection) {
		CD.resize(nodes->size(), 4);
		CDe.resize(2, 2);
	}

	const MaterialBaseConfiguration *phase1, *phase2;
	if (material->phase_change) {
		phase1 = &material->phases.find(1)->second;
		phase2 = &material->phases.find(2)->second;
	}

	for (size_t n = 0; n < nodes->size(); n++) {
		auto it = std::lower_bound(intervals.begin(), intervals.end(), nodes->at(n), [] (const DomainInterval &interval, eslocal node) { return interval.end < node; });
		T(n, 0) = (*_temperature->decomposedData)[domain][it->DOFOffset + nodes->at(n) - it->begin];
		const Point &p = _mesh->nodes->coordinates->datatarray()[nodes->at(n)];
		coordinates(n, 0) = p.x;
		coordinates(n, 1) = p.y;
		center += p;
		thickness(n, 0) = thick != NULL ? thick->evaluate(p, T(n, 0), _step->currentTime) : 1;
		if (material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, material->phase_change_temperature - material->transition_interval / 2, material->phase_change_temperature + material->transition_interval / 2, T(n, 0), material->smooth_step_order);
			assembleMaterialMatrix(n, p, phase1, phase, T(n, 0), K, CD, tangentCorrection);
			assembleMaterialMatrix(n, p, phase2, (1 - phase), T(n, 0), K, CD, tangentCorrection);
			m(n, 0) =
					(    phase  * phase1->density.evaluator->evaluate(p, _step->currentTime, T(n, 0)) +
					(1 - phase) * phase2->density.evaluator->evaluate(p, _step->currentTime, T(n, 0))) *

					(    phase  * phase1->heat_capacity.evaluator->evaluate(p, _step->currentTime, T(n, 0)) +
					(1 - phase) * phase2->heat_capacity.evaluator->evaluate(p, _step->currentTime, T(n, 0)) +
					material->latent_heat * derivation) * thickness(n, 0);
		} else {
			assembleMaterialMatrix(n, p, material, 1, T(n, 0), K, CD, tangentCorrection);
			m(n, 0) =
					material->density.evaluator->evaluate(p, _step->currentTime, T(n, 0)) *
					material->heat_capacity.evaluator->evaluate(p, _step->currentTime, T(n, 0)) * thickness(n, 0);
		}

		if (translation_motion) {
			U(n, 0) = translation_motion->x.evaluator->evaluate(p, _step->currentTime, T(n, 0)) * m(n, 0);
			U(n, 1) = translation_motion->y.evaluator->evaluate(p, _step->currentTime, T(n, 0)) * m(n, 0);
		}
		if (heat_source) {
			f(n, 0) = heat_source->evaluate(p, _step->currentTime, T(n, 0)) * thickness(n, 0);
		}
	}

	center /= nodes->size();
	if (translation_motion) {
		(*_translationMotion->data)[2 * eindex + 0] = translation_motion->x.evaluator->evaluate(center, _step->currentTime, 0);
		(*_translationMotion->data)[2 * eindex + 1] = translation_motion->y.evaluator->evaluate(center, _step->currentTime, 0);
	}

	eslocal Ksize = nodes->size();

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if ((matrices & Matrices::K) || ((matrices & Matrices::R) && _step->timeIntegrationConstantK != 0)) {
		Ke.resize(Ksize, Ksize);
		Ke = 0;
	}
	if ((matrices & Matrices::M) || ((matrices & Matrices::R) && _step->timeIntegrationConstantM != 0)) {
		Me.resize(Ksize, Ksize);
		Me = 0;
	}
	if (matrices & Matrices::R) {
		Re.resize(Ksize, 1);
		Re = 0;
	}
	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	if (tangentCorrection) {
		tangentK.resize(Ksize, Ksize);
	}

	DenseMatrix g(1, 2), u(1, 2), v(1, 2), re(1, nodes->size());
	double normGradN = 0;

	if ((matrices & Matrices::M) && _configuration.diffusion_split) {
		g(0, 0) = (*_gradient->data)[2 * eindex + 0];
		g(0, 1) = (*_gradient->data)[2 * eindex + 1];
	}

//	std::cout << "C: " << coordinates;
//	std::cout << "K: " << K;

//	if (eindex == 0) {
//		std::cout << "C: " << coordinates;
//		std::cout << "K: " << K;
//	}

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

		DenseMatrix b_e(1, nodes->size()), b_e_c(1, nodes->size()), g_e(1, nodes->size());
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

		if (_configuration.diffusion_split && g.norm() != 0) {
			gh_e = 2 * g.norm() / g_e.norm();
			tauK = (C1 * gh_e * gh_e) / (Ce(0, 0) * C2 + gh_e * gh_e * (gpM(0, 0) / _step->timeStep));
			xi = std::max(1., 1 / (1 - tauK * gpM(0, 0) / _step->timeStep));
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

		Ce(0, 0) += _configuration.sigma * h_e * norm_u_e;
		Ce(1, 1) += _configuration.sigma * h_e * norm_u_e;

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
			if (_configuration.diffusion_split) {
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
			for (eslocal i = 0; i < Ksize; i++) {
				fe(i, 0) += detJ * weighFactor[gp] * N[gp](0, i) * f(i, 0);
				if (norm_u_e != 0) {
					fe(i, 0) += detJ * weighFactor[gp] * h_e * tau_e * b_e(0, i) * f(i, 0) / (2 * norm_u_e);
				}
			}
		}
	}

	if ((matrices & Matrices::M) && _configuration.diffusion_split) {
		DenseMatrix T1, T2;
		T1.multiply(Ke, T, 1, 0);
		T2.multiply(gKe, T, 1, 0);
		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += T1(i, 0) - T2(i, 0);
		}
	}

	if (matrices & Matrices::R) {
		Re.multiply(Ke, T, _step->timeIntegrationConstantK, 0);
		Re.multiply(Me, T, _step->timeIntegrationConstantM, 1);
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

//	if (eindex == 0) {
//		std::cout << " -- " << eindex << " -- \n" << Ke << fe;
//	}
}

void HeatTransfer2D::processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	ESINFO(ERROR) << "Advection diffusion 2D cannot process face";
}


void HeatTransfer2D::processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	const ConvectionConfiguration *convection = NULL;
	const RadiationConfiguration *radiation = NULL;
	const Evaluator *heatFlow = NULL, *heatFlux = NULL, *thick = NULL;

	auto itc = _configuration.load_steps_settings.at(_step->step + 1).convection.find(region->name);
	if (itc != _configuration.load_steps_settings.at(_step->step + 1).convection.end()) {
		convection = &itc->second;
	}

	auto itr = _configuration.load_steps_settings.at(_step->step + 1).diffuse_radiation.find(region->name);
	if (itr != _configuration.load_steps_settings.at(_step->step + 1).diffuse_radiation.end()) {
		radiation = &itr->second;
	}

	auto it = _configuration.load_steps_settings.at(_step->step + 1).heat_flow.find(region->name);
	if (it != _configuration.load_steps_settings.at(_step->step + 1).heat_flow.end()) {
		heatFlow = it->second.evaluator;
	}
	it = _configuration.load_steps_settings.at(_step->step + 1).heat_flux.find(region->name);
	if (it != _configuration.load_steps_settings.at(_step->step + 1).heat_flux.end()) {
		heatFlux = it->second.evaluator;
	}

	for (auto it = _configuration.thickness.begin(); it != _configuration.thickness.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			thick = it->second.evaluator;
			break;
		}
	}

	if (convection == NULL && heatFlow == NULL && heatFlux == NULL && radiation == NULL) {
		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}
	if (!(matrices & (Matrices::K | Matrices::f))) {
		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}

	auto nodes = region->elements->cbegin() + eindex;
	auto epointer = region->epointers->datatarray()[eindex];
	const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[domain];

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	DenseMatrix coordinates(nodes->size(), 2), dND(1, 2), q(nodes->size(), 1), htc(nodes->size(), 1), thickness(nodes->size(), 1), flow(nodes->size(), 1), emiss(nodes->size(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpThickness(1, 1), gpFlow(1, 1), gpEmiss(1, 1);

	double area = region->area, temp, text;
	eslocal Ksize = nodes->size();
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);

	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	if (convection != NULL || radiation != NULL) {
		Ke.resize(Ksize, Ksize);
		Ke = 0;
	}

	for (size_t n = 0; n < nodes->size(); n++) {
		auto it = std::lower_bound(intervals.begin(), intervals.end(), nodes->at(n), [] (const DomainInterval &interval, eslocal node) { return interval.end < node; });
		temp = (*_temperature->decomposedData)[domain][it->DOFOffset + nodes->at(n) - it->begin];
		const Point &p = _mesh->nodes->coordinates->datatarray()[nodes->at(n)];
		coordinates(n, 0) = p.x;
		coordinates(n, 1) = p.y;

		if (convection != NULL) {
			text = convection->external_temperature.evaluator->evaluate(p, temp, _step->currentTime);
			htc(n, 0) = computeHTC(convection, p, temp);

			if (_step->iteration) {
				q(n, 0) += htc(n, 0) * (text - temp);
			} else {
				q(n, 0) += htc(n, 0) * (text);
			}
		}

		if (radiation != NULL) {
			emiss(n, 0) = CONST_Stefan_Boltzmann * radiation->emissivity.evaluator->evaluate(p, temp, _step->currentTime);
			q(n, 0) += emiss(n, 0) * (pow(radiation->external_temperature.evaluator->evaluate(p, temp, _step->currentTime), 4) - pow(temp, 4));
			emiss(n, 0) *= 4 * temp * temp * temp;
		}
		if (heatFlow) {
			q(n, 0) += heatFlow->evaluate(p, temp, _step->currentTime) / area;
		}
		if (heatFlux) {
			q(n, 0) += heatFlux->evaluate(p, temp, _step->currentTime);
		}

		thickness(n, 0) = thick != NULL ? thick->evaluate(p, temp, _step->currentTime) : 1;
		q(n, 0) *= thickness(n, 0);
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);

		if (convection != NULL || radiation != NULL) {
			gpThickness.multiply(N[gp], thickness);
		}
		if (convection != NULL) {
			gpHtc.multiply(N[gp], htc);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
		}
		if (radiation != NULL) {
			gpEmiss.multiply(N[gp], emiss);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0), 1, true);
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % nodes->size()) * gpQ(0, 0);
		}
	}
}

void HeatTransfer2D::processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
}

void HeatTransfer2D::postProcessElement(eslocal domain, eslocal eindex)
{
	auto nodes = _mesh->elements->procNodes->cbegin() + eindex;
	auto epointer = _mesh->elements->epointers->datatarray()[eindex];
	const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[domain];
	const ECFExpressionVector *translation_motion = NULL;
//	Evaluator *heat_source = NULL;
	Evaluator *thick = NULL;
	for (auto it = _configuration.load_steps_settings.at(_step->step + 1).translation_motions.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).translation_motions.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			translation_motion = &it->second;
			break;
		}
	}
//	for (auto it = _configuration.load_steps_settings.at(_step->step + 1).heat_source.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).heat_source.end(); ++it) {
//		ElementsRegionStore *region = _mesh->eregion(it->first);
//		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
//			heat_source = it->second.evaluator;
//			break;
//		}
//	}
	for (auto it = _configuration.thickness.begin(); it != _configuration.thickness.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			thick = it->second.evaluator;
			break;
		}
	}

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
//	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	DenseMatrix Ce(2, 2), coordinates, J(2, 2), invJ(2, 2), dND, temp(nodes->size(), 1);
	double detJ, m, norm_u_e, h_e;
	DenseMatrix thickness(nodes->size(), 1), U(nodes->size(), 2), K(nodes->size(), 4), gpK(1, 4), CD;
	DenseMatrix u(1, 2), matFlux(2, 1), matGradient(2, 1);

	const MaterialConfiguration* material = _mesh->materials[_mesh->elements->material->datatarray()[eindex]];

	const MaterialBaseConfiguration *phase1, *phase2;
	if (material->phase_change) {
		phase1 = &material->phases.find(1)->second;
		phase2 = &material->phases.find(2)->second;
	}

	coordinates.resize(nodes->size(), 2);

	for (size_t n = 0; n < nodes->size(); n++) {
		auto it = std::lower_bound(intervals.begin(), intervals.end(), nodes->at(n), [] (const DomainInterval &interval, eslocal node) { return interval.end < node; });
		temp(n, 0) = (*_temperature->decomposedData)[domain][it->DOFOffset + nodes->at(n) - it->begin];
		const Point &p = _mesh->nodes->coordinates->datatarray()[nodes->at(n)];
		coordinates(n, 0) = p.x;
		coordinates(n, 1) = p.y;
		thickness(n, 0) = thick != NULL ? thick->evaluate(p, temp(n, 0), _step->currentTime) : 1;
		if (material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, material->phase_change_temperature - material->transition_interval / 2, material->phase_change_temperature + material->transition_interval / 2, temp(n, 0), material->smooth_step_order);
			assembleMaterialMatrix(n, p, phase1, phase, temp(n, 0), K, CD, false);
			assembleMaterialMatrix(n, p, phase2, (1 - phase), temp(n, 0), K, CD, false);
			m =
					(    phase  * phase1->density.evaluator->evaluate(p, _step->currentTime, temp(n, 0)) +
					(1 - phase) * phase2->density.evaluator->evaluate(p, _step->currentTime, temp(n, 0))) *

					(    phase  * phase1->heat_capacity.evaluator->evaluate(p, _step->currentTime, temp(n, 0)) +
					(1 - phase) * phase2->heat_capacity.evaluator->evaluate(p, _step->currentTime, temp(n, 0)) +
					material->latent_heat * derivation) * thickness(n, 0);
			if (_phaseChange) {
				(*_phaseChange->decomposedData)[domain][it->DOFOffset + nodes->at(n) - it->begin] = phase;
				(*_latentHeat->decomposedData)[domain][it->DOFOffset + nodes->at(n) - it->begin] = material->latent_heat * derivation;
			}
		} else {
			assembleMaterialMatrix(n, p, material, 1, temp(n, 0), K, CD, false);
			m =
					material->density.evaluator->evaluate(p, _step->currentTime, temp(n, 0)) *
					material->heat_capacity.evaluator->evaluate(p, _step->currentTime, temp(n, 0)) * thickness(n, 0);
		}

		if (translation_motion) {
			U(n, 0) = translation_motion->x.evaluator->evaluate(p, _step->currentTime, temp(n, 0)) * m;
			U(n, 1) = translation_motion->y.evaluator->evaluate(p, _step->currentTime, temp(n, 0)) * m;
		}
//		if (heat_source) {
//			f(n, 0) = heat_source->evaluate(p, _step->currentTime, temp(n, 0));
//		}
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
			DenseMatrix b_e(1, nodes->size());
			b_e.multiply(u, dND, 1, 0);
			h_e = 2 * norm_u_e / b_e.norm();
		}

		Ce(0, 0) += _configuration.sigma * h_e * norm_u_e;
		Ce(1, 1) += _configuration.sigma * h_e * norm_u_e;

		if (_propertiesConfiguration.gradient) {
			matGradient.multiply(dND, temp, 1, 1);
		}
		if (_propertiesConfiguration.flux) {
			matFlux.multiply(Ce, dND * temp, 1, 1);
		}
	}

	if (_propertiesConfiguration.gradient) {
		(*_gradient->data)[2 * eindex + 0] = matGradient(0, 0) / N.size();
		(*_gradient->data)[2 * eindex + 1] = matGradient(1, 0) / N.size();
	}

	if (_propertiesConfiguration.flux) {
		(*_flux->data)[2 * eindex + 0] = matFlux(0, 0) / N.size();
		(*_flux->data)[2 * eindex + 1] = matFlux(1, 0) / N.size();
	}
}

void HeatTransfer2D::processSolution()
{
	if (_gradient || _flux || _phaseChange) {
		#pragma omp parallel for
		for (eslocal d = 0; d < _mesh->elements->ndomains; d++) {
			for (eslocal e = _mesh->elements->elementsDistribution[d]; e < (eslocal)_mesh->elements->elementsDistribution[d + 1]; e++) {
				postProcessElement(d, e);
			}
		}
	}
}

