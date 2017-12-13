
#include "../../config/ecf/physics/heattransfer.h"
#include "../../config/ecf/output.h"
#include "../step.h"
#include "../instance.h"
#include "../solution.h"
#include "../constraints/equalityconstraints.h"

#include "../../mesh/mesh.h"


#include "../../basis/matrices/denseMatrix.h"
#include "../../solver/generic/SparseMatrix.h"
#include "heattransfer2d.h"

using namespace espreso;

HeatTransfer2D::HeatTransfer2D(Mesh *mesh, Instance *instance, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration)
: Physics("HEAT TRANSFER 2D", mesh, instance, &configuration), HeatTransfer(configuration, propertiesConfiguration)
{
	_equalityConstraints = new EqualityConstraints(*_instance, *_mesh, {}, 1, configuration.load_steps_settings.at(1).feti.redundant_lagrange, configuration.load_steps_settings.at(1).feti.scaling);
}

void HeatTransfer2D::assembleMaterialMatrix(const Step &step, eslocal eindex, eslocal node, const MaterialBaseConfiguration *mat, double phase, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const
{
//	auto conductivity = [&] (int row, int column, double t) {
//		return mat->thermal_conductivity.values.get(row, column).evaluate(_mesh->coordinates()[e->node(node)], step.currentTime, t);
//	};
//
//	auto derivation = [&] (int row, int column, double h) {
//		return (conductivity(row, column, temp + h) - conductivity(row, column, temp - h)) / (2 * h);
//	};
//
//	auto d2r = [] (double degree) -> double {
//		return M_PI * degree / 180;
//	};
//
//	double cos, sin;
//	switch (mat->coordinate_system.type) {
//	case CoordinateSystemConfiguration::TYPE::CARTESIAN:
//		cos = std::cos(d2r(mat->coordinate_system.rotation_z.evaluate(_mesh->coordinates()[e->node(node)], step.currentTime, temp)));
//		sin = std::sin(d2r(mat->coordinate_system.rotation_z.evaluate(_mesh->coordinates()[e->node(node)], step.currentTime, temp)));
//		break;
//	case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: {
//		Point origin(
//				mat->coordinate_system.center_x.evaluate(_mesh->coordinates()[e->node(node)], step.currentTime, temp),
//				mat->coordinate_system.center_y.evaluate(_mesh->coordinates()[e->node(node)], step.currentTime, temp),
//				0);
//		const Point &p = _mesh->coordinates()[e->node(node)];
//		double rotation = std::atan2((p.y - origin.y), (p.x - origin.x));
//		cos = std::cos(rotation);
//		sin = std::sin(rotation);
//		break;
//	}
//	default:
//		ESINFO(ERROR) << "Invalid material type (SPHERICAL for 2D).";
//	}
//
//	DenseMatrix TCT(2, 2), T(2, 2), C(2, 2), _CD, TCDT;
//	T(0, 0) =  cos; T(0, 1) = sin;
//	T(1, 0) = -sin; T(1, 1) = cos;
//
//	if (tangentCorrection) {
//		_CD.resize(2, 2);
//		TCDT.resize(2, 2);
//	}
//
//
//
//	switch (mat->thermal_conductivity.model) {
//	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
//		C(0, 0) = C(1, 1) = conductivity(0, 0, temp);
//		C(0, 1) = C(1, 0) = 0;
//		if (tangentCorrection) {
//			_CD(0, 0) = _CD(1, 1) = derivation(0, 0, temp / 1e4);
//			_CD(0, 1) = _CD(1, 0) = 0;
//		}
//		break;
//	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
//		C(0, 0) = conductivity(0, 0, temp);
//		C(1, 1) = conductivity(1, 1, temp);
//		C(0, 1) = C(1, 0) = 0;
//		if (tangentCorrection) {
//			_CD(0, 0) = derivation(0, 0, temp / 1e4);
//			_CD(1, 1) = derivation(1, 1, temp / 1e4);
//			_CD(0, 1) = _CD(1, 0) = 0;
//		}
//		break;
//	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
//		C(0, 0) = conductivity(0, 0, temp);
//		C(1, 1) = conductivity(1, 1, temp);
//		C(1, 0) = C(0, 1) = conductivity(0, 1, temp);
//		if (tangentCorrection) {
//			_CD(0, 0) = derivation(0, 0, temp / 1e4);
//			_CD(1, 1) = derivation(1, 1, temp / 1e4);
//			_CD(0, 1) = _CD(1, 0) = derivation(0, 1, temp / 1e4);
//		}
//		break;
//	case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
//		C(0, 0) = conductivity(0, 0, temp);
//		C(1, 1) = conductivity(1, 1, temp);
//		C(0, 1) = conductivity(0, 1, temp);
//		C(1, 0) = conductivity(1, 0, temp);
//		if (tangentCorrection) {
//			_CD(0, 0) = derivation(0, 0, temp / 1e4);
//			_CD(1, 1) = derivation(1, 1, temp / 1e4);
//			_CD(0, 1) = derivation(0, 1, temp / 1e4);
//			_CD(1, 0) = derivation(1, 0, temp / 1e4);
//		}
//		break;
//	default:
//		ESINFO(ERROR) << "Advection diffusion 2D not supports set material model";
//	}
//
//	if (tangentCorrection) {
//		TCDT.multiply(T, _CD * T, 1, 0, true, false);
//		CD(node, 0) += phase * TCDT(0, 0);
//		CD(node, 1) += phase * TCDT(1, 1);
//		CD(node, 2) += phase * TCDT(0, 1);
//		CD(node, 3) += phase * TCDT(1, 0);
//	}
//
//	TCT.multiply(T, C * T, 1, 0, true, false);
//	K(node, 0) += phase * TCT(0, 0);
//	K(node, 1) += phase * TCT(1, 1);
//	K(node, 2) += phase * TCT(0, 1);
//	K(node, 3) += phase * TCT(1, 0);
}

void HeatTransfer2D::processElement(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
//	bool CAU = _configuration.stabilization == HeatTransferConfiguration::STABILIZATION::CAU;
//	bool tangentCorrection = (matrices & Matrices::K) && step.tangentMatrixCorrection;
//
//	DenseMatrix Ce(2, 2), coordinates(e->nodes(), 2), J(2, 2), invJ(2, 2), dND;
//	double detJ, temp;
//	DenseMatrix f(e->nodes(), 1);
//	DenseMatrix U(e->nodes(), 2);
//	DenseMatrix m(e->nodes(), 1);
//	DenseMatrix T(e->nodes(), 1);
//	DenseMatrix thickness(e->nodes(), 1), K(e->nodes(), 4);
//	DenseMatrix gpThickness(1, 1), gpK(1, 4), gpM(1, 1);
//	DenseMatrix tangentK, BT, BTN, gpCD, CD, CDBTN, CDe;

//	const MaterialConfiguration* material = _mesh->materials()[e->param(OldElement::MATERIAL)];
//
//	if (tangentCorrection) {
//		CD.resize(e->nodes(), 4);
//		CDe.resize(2, 2);
//	}
//
//	const MaterialBaseConfiguration *phase1, *phase2;
//	if (material->phase_change) {
//		phase1 = &material->phases.find(1)->second;
//		phase2 = &material->phases.find(2)->second;
//	}
//
//	for (size_t i = 0; i < e->nodes(); i++) {
//		temp = solution[offset + SolutionIndex::TEMPERATURE]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
//		T(i, 0) = temp;
//		coordinates(i, 0) = _mesh->coordinates()[e->node(i)].x;
//		coordinates(i, 1) = _mesh->coordinates()[e->node(i)].y;
//		thickness(i, 0) = e->getProperty(Property::THICKNESS, step.step, _mesh->coordinates()[e->node(i)], step.currentTime, temp, 1);
//		if (material->phase_change) {
//			double phase, derivation;
//			smoothstep(phase, derivation, material->phase_change_temperature - material->transition_interval / 2, material->phase_change_temperature + material->transition_interval / 2, temp, material->smooth_step_order);
//			assembleMaterialMatrix(step, e, i, phase1, phase, temp, K, CD, tangentCorrection);
//			assembleMaterialMatrix(step, e, i, phase2, (1 - phase), temp, K, CD, tangentCorrection);
//			m(i, 0) =
//					(    phase  * phase1->density.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp) +
//					(1 - phase) * phase2->density.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp)) *
//
//					(    phase  * phase1->heat_capacity.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp) +
//					(1 - phase) * phase2->heat_capacity.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp) +
//					material->latent_heat * derivation) *
//
//					thickness(i, 0);
//		} else {
//			assembleMaterialMatrix(step, e, i, material, 1, temp, K, CD, tangentCorrection);
//			m(i, 0) =
//					material->density.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp) *
//					material->heat_capacity.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp) *
//					thickness(i, 0);
//		}
//
//		U(i, 0) = e->getProperty(Property::TRANSLATION_MOTION_X, step.step, _mesh->coordinates()[e->node(i)], step.currentTime, temp, 0) * m(i, 0);
//		U(i, 1) = e->getProperty(Property::TRANSLATION_MOTION_Y, step.step, _mesh->coordinates()[e->node(i)], step.currentTime, temp, 0) * m(i, 0);
//		f(i, 0) = e->sumProperty(Property::HEAT_SOURCE, step.step, _mesh->coordinates()[e->node(i)], step.currentTime, temp, 0) * thickness(i, 0);
//	}
//
//	eslocal Ksize = e->nodes();
//
//	Ke.resize(0, 0);
//	Me.resize(0, 0);
//	Re.resize(0, 0);
//	fe.resize(0, 0);
//	if ((matrices & Matrices::K) || ((matrices & Matrices::R) && step.timeIntegrationConstantK != 0)) {
//		Ke.resize(Ksize, Ksize);
//		Ke = 0;
//	}
//	if ((matrices & Matrices::M) || ((matrices & Matrices::R) && step.timeIntegrationConstantM != 0)) {
//		Me.resize(Ksize, Ksize);
//		Me = 0;
//	}
//	if (matrices & Matrices::R) {
//		Re.resize(Ksize, 1);
//		Re = 0;
//	}
//	if (matrices & Matrices::f) {
//		fe.resize(Ksize, 1);
//		fe = 0;
//	}
//
//	if (tangentCorrection) {
//		tangentK.resize(Ksize, Ksize);
//	}
//
//	DenseMatrix u(1, 2), v(1, 2), re(1, e->nodes());
//	double normGradN = 0;
//
//	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
//		u.multiply(e->N()[gp], U, 1, 0);
//
//		J.multiply(e->dN()[gp], coordinates);
//		detJ = determinant2x2(J.values());
//		inverse2x2(J.values(), invJ.values(), detJ);
//
//		gpThickness.multiply(e->N()[gp], thickness);
//		gpK.multiply(e->N()[gp], K);
//		if (tangentCorrection) {
//			gpCD.multiply(e->N()[gp], CD);
//			CDe(0, 0) = gpCD(0, 0);
//			CDe(1, 1) = gpCD(0, 1);
//			CDe(0, 1) = gpCD(0, 2);
//			CDe(1, 0) = gpCD(0, 3);
//		}
//		gpM.multiply(e->N()[gp], m);
//
//		Ce(0, 0) = gpK(0, 0);
//		Ce(1, 1) = gpK(0, 1);
//		Ce(0, 1) = gpK(0, 2);
//		Ce(1, 0) = gpK(0, 3);
//
//
//		dND.multiply(invJ, e->dN()[gp]);
//
//		DenseMatrix b_e(1, e->nodes()), b_e_c(1, e->nodes());
//		b_e.multiply(u, dND, 1, 0);
//
//		if (CAU) {
//			normGradN = dND.norm();
//			if (normGradN >= 1e-12) {
//				for (size_t i = 0; i < re.columns(); i++) {
//					re(0, i) = b_e(0, i) - f(i, 0);
//				}
//				DenseMatrix ReBt(1, 2);
//				ReBt.multiply(re, dND, 1 / pow(normGradN, 2), 0, false, true);
//				for (size_t i = 0; i < ReBt.columns(); i++) {
//					v(0, i) = u(0, i) - ReBt(0, i);
//				}
//			} else {
//				v = u;
//			}
//		}
//
//
//		double norm_u_e = u.norm();
//		double h_e = 0, tau_e = 0, konst = 0;
//		double C_e = 0;
//
//		if (norm_u_e != 0) {
//			h_e = 2 * norm_u_e / b_e.norm();
//			double P_e = h_e * norm_u_e / (2 * Ce(0, 0));
//			tau_e = std::max(0.0, 1 - 1 / P_e);
//			konst = h_e * tau_e / (2 * norm_u_e);
//
//			if (CAU) {
//				DenseMatrix u_v(1, 2);
//				u_v(0, 0) = u(0, 0) - v(0, 0);
//				u_v(0, 1) = u(0, 1) - v(0, 1);
//				b_e_c.multiply(u_v, dND, 1, 0);
//				double norm_u_v = u_v.norm();
//				double h_e_c = 2 * norm_u_v / b_e_c.norm();
//				double P_e_c = h_e_c * norm_u_v / (2 * Ce.norm());
//				double tau_e_c = std::max(0.0, 1 - 1 / P_e_c);
//
//				double konst1 = re.norm() / normGradN;
//				double konst2 = tau_e * h_e != 0 ? tau_e_c * h_e_c / (tau_e * h_e) : 0;
//				if (konst1 / norm_u_e < konst2) {
//					C_e = tau_e * h_e * konst1 * (konst2 - konst1 / norm_u_e) / 2;
//				} else {
//					C_e = 0;
//				}
//			}
//		}
//
//		Ce(0, 0) += _configuration.sigma * h_e * norm_u_e;
//		Ce(1, 1) += _configuration.sigma * h_e * norm_u_e;
//
//		if (matrices & (Matrices::M | Matrices::R)) {
//			Me.multiply(e->N()[gp], e->N()[gp], detJ * gpM(0, 0) * e->weighFactor()[gp], 1, true);
//		}
//		if (matrices & (Matrices::K | Matrices::R)) {
//			if (tangentCorrection) {
//				BT.multiply(dND, T);
//				BTN.multiply(BT, e->N()[gp]);
//				CDBTN.multiply(CDe, BTN);
//				tangentK.multiply(dND, CDBTN,  detJ * e->weighFactor()[gp] * gpThickness(0, 0), 1, true);
//			}
//			Ke.multiply(dND, Ce * dND, detJ * e->weighFactor()[gp] * gpThickness(0, 0), 1, true);
//			Ke.multiply(e->N()[gp], b_e, detJ * e->weighFactor()[gp], 1, true);
//			if (konst * e->weighFactor()[gp] * detJ != 0) {
//				Ke.multiply(b_e, b_e, konst * e->weighFactor()[gp] * detJ, 1, true);
//			}
//			if (CAU) {
//				Ke.multiply(dND, dND, C_e * e->weighFactor()[gp] * detJ, 1, true);
//			}
//		}
//
//		if (matrices & Matrices::f) {
//			for (eslocal i = 0; i < Ksize; i++) {
//				fe(i, 0) += detJ * e->weighFactor()[gp] * e->N()[gp](0, i) * f(i, 0);
//				if (norm_u_e != 0) {
//					fe(i, 0) += detJ * e->weighFactor()[gp] * h_e * tau_e * b_e(0, i) * f(i, 0) / (2 * norm_u_e);
//				}
//			}
//		}
//	}
//
//	if (matrices & Matrices::R) {
//		Re.multiply(Ke, T, step.timeIntegrationConstantK, 0);
//		Re.multiply(Me, T, step.timeIntegrationConstantM, 1);
//		if (!(matrices & Matrices::K)) {
//			Ke.resize(0, 0);
//		}
//		if (!(matrices & Matrices::M)) {
//			Me.resize(0, 0);
//		}
//	}
//
//	if (tangentCorrection) {
//		Ke += tangentK;
//	}
}

void HeatTransfer2D::processFace(const Step &step, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	ESINFO(ERROR) << "Advection diffusion 2D cannot process face";
}


void HeatTransfer2D::processEdge(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
//	if (!(e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step) ||
//		e->hasProperty(Property::HEAT_FLOW, step.step) ||
//		e->hasProperty(Property::HEAT_FLUX, step.step))) {
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
//	DenseMatrix coordinates(e->nodes(), 2), dND(1, 2), q(e->nodes(), 1), htc(e->nodes(), 1), thickness(e->nodes(), 1), flow(e->nodes(), 1), emiss(e->nodes(), 1);
//	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpThickness(1, 1), gpFlow(1, 1), gpEmiss(1, 1);
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
//		if (step.step < e->regions()[r]->settings.size() && e->regions()[r]->settings[step.step].count(Property::HEAT_FLOW)) {
//			area = e->regions()[r]->area;
//			break;
//		}
//	}
//	if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step)) {
//		Ke.resize(Ksize, Ksize);
//		Ke = 0;
//	}
//
//	const std::vector<DenseMatrix> &dN = e->dN();
//	const std::vector<DenseMatrix> &N = e->N();
//	const std::vector<double> &weighFactor = e->weighFactor();
//
//
//	const ConvectionConfiguration *convection = NULL;
//	for (size_t r = 0; convection == NULL && r < e->regions().size(); r++) {
//		auto regionit = _configuration.load_steps_settings.at(step.step + 1).convection.find(e->regions()[r]->name);
//		if (regionit != _configuration.load_steps_settings.at(step.step + 1).convection.end()) {
//			convection = &regionit->second;
//		}
//	}

//	for (size_t n = 0; n < e->nodes(); n++) {
//		coordinates(n, 0) = _mesh->coordinates()[e->node(n)].x;
//		coordinates(n, 1) = _mesh->coordinates()[e->node(n)].y;
//
//		temp = solution[offset + SolutionIndex::TEMPERATURE]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(n), e->domains().front()));
//		htc(n, 0) = convection != NULL ? computeHTC(*convection, e, _mesh->coordinates()[e->node(n)], step, temp) : 0;
//
//
//		if (step.iteration) {
//			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, step.step, _mesh->coordinates()[e->node(n)], step.currentTime, temp, 0) - temp);
//		} else {
//			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, step.step, _mesh->coordinates()[e->node(n)], step.currentTime, temp, 0));
//		}
//
//		emiss(n, 0) = CONST_Stefan_Boltzmann * e->getProperty(Property::EMISSIVITY, step.step, _mesh->coordinates()[e->node(n)], step.currentTime, temp, 0);
//		q(n, 0) += emiss(n, 0) * (pow(e->getProperty(Property::EXTERNAL_TEMPERATURE, step.step, _mesh->coordinates()[e->node(n)], step.currentTime, temp, 0), 4) - pow(temp, 4));
//		q(n, 0) += e->getProperty(Property::HEAT_FLOW, step.step, _mesh->coordinates()[e->node(n)], step.currentTime, temp, 0) / area;
//		q(n, 0) += e->getProperty(Property::HEAT_FLUX, step.step, _mesh->coordinates()[e->node(n)], step.currentTime, temp, 0);
//
//		emiss(n, 0) *= 4 * temp * temp * temp;
//
//		thickness(n, 0) = e->getProperty(Property::THICKNESS, step.step, _mesh->coordinates()[e->node(n)], step.currentTime, temp, 1);
//		q(n, 0) *= thickness(n, 0);
//	}
//
//	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
//		dND.multiply(dN[gp], coordinates);
//		double J = dND.norm();
//		gpQ.multiply(N[gp], q);
//		if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step)) {
//			gpHtc.multiply(N[gp], htc);
//			gpEmiss.multiply(N[gp], emiss);
//			gpThickness.multiply(N[gp], thickness);
//
//			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0) * gpThickness(0, 0), 1, true);
//			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0) * gpThickness(0, 0), 1, true);
//		}
//		for (eslocal i = 0; i < Ksize; i++) {
//			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % e->nodes()) * gpQ(0, 0);
//		}
//	}
}

void HeatTransfer2D::processNode(const Step &step, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
}

void HeatTransfer2D::postProcessElement(const Step &step, eslocal eindex, std::vector<Solution*> &solution)
{
//	DenseMatrix Ce(2, 2), coordinates, J(2, 2), invJ(2, 2), dND, temp(e->nodes(), 1);
//	double detJ, m, norm_u_e, h_e;
//	DenseMatrix thickness(e->nodes(), 1), U(e->nodes(), 2), K(e->nodes(), 4), gpK(1, 4), CD;
//	DenseMatrix u(1, 2), matFlux(2, 1), matGradient(2, 1);

//	const MaterialConfiguration* material = _mesh->materials()[e->param(OldElement::MATERIAL)];
//
//	const MaterialBaseConfiguration *phase1, *phase2;
//	if (material->phase_change) {
//		phase1 = &material->phases.find(1)->second;
//		phase2 = &material->phases.find(2)->second;
//	}
//
//	coordinates.resize(e->nodes(), 2);
//
//	for (size_t i = 0; i < e->nodes(); i++) {
//		temp(i, 0) = solution[offset + SolutionIndex::TEMPERATURE]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
//		coordinates(i, 0) = _mesh->coordinates()[e->node(i)].x;
//		coordinates(i, 1) = _mesh->coordinates()[e->node(i)].y;
//		thickness(i, 0) = e->getProperty(Property::THICKNESS, step.step, _mesh->coordinates()[e->node(i)], step.currentTime, temp(i, 0), 1);
//		if (material->phase_change) {
//			double phase, derivation;
//			smoothstep(phase, derivation, material->phase_change_temperature - material->transition_interval / 2, material->phase_change_temperature + material->transition_interval / 2, temp(i, 0), material->smooth_step_order);
//			assembleMaterialMatrix(step, e, i, phase1, phase, temp(i, 0), K, CD, false);
//			assembleMaterialMatrix(step, e, i, phase2, (1 - phase), temp(i, 0), K, CD, false);
//			m =
//					(    phase  * phase1->density.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp(i, 0)) +
//					(1 - phase) * phase2->density.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp(i, 0))) *
//
//					(    phase  * phase1->heat_capacity.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp(i, 0)) +
//					(1 - phase) * phase2->heat_capacity.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp(i, 0)) +
//					material->latent_heat * derivation) *
//
//					thickness(i, 0);
//			solution[offset + SolutionIndex::PHASE]->data[e->domains().front()][_mesh->coordinates().localIndex(e->node(i), e->domains().front())] = phase;
//			solution[offset + SolutionIndex::LATENT_HEAT]->data[e->domains().front()][_mesh->coordinates().localIndex(e->node(i), e->domains().front())] = material->latent_heat * derivation;
//		} else {
//			assembleMaterialMatrix(step, e, i, material, 1, temp(i, 0), K, CD, false);
//			m =
//					material->density.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp(i, 0)) *
//					material->heat_capacity.evaluate(_mesh->coordinates()[e->node(i)], step.currentTime, temp(i, 0)) *
//					thickness(i, 0);
//		}
//
//		U(i, 0) = e->getProperty(Property::TRANSLATION_MOTION_X, step.step, _mesh->coordinates()[e->node(i)], step.currentTime, temp(i, 0), 0) * m;
//		U(i, 1) = e->getProperty(Property::TRANSLATION_MOTION_Y, step.step, _mesh->coordinates()[e->node(i)], step.currentTime, temp(i, 0), 0) * m;
//	}
//
//
//
//	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
//		u.multiply(e->N()[gp], U, 1, 0);
//
//		J.multiply(e->dN()[gp], coordinates);
//		detJ = determinant2x2(J.values());
//		inverse2x2(J.values(), invJ.values(), detJ);
//
//		gpK.multiply(e->N()[gp], K);
//
//		Ce(0, 0) = gpK(0, 0);
//		Ce(1, 1) = gpK(0, 1);
//		Ce(0, 1) = gpK(0, 2);
//		Ce(1, 0) = gpK(0, 3);
//
//		dND.multiply(invJ, e->dN()[gp]);
//
//		norm_u_e = u.norm();
//		h_e = 0;
//
//		if (norm_u_e != 0) {
//			DenseMatrix b_e(1, e->nodes());
//			b_e.multiply(u, dND, 1, 0);
//			h_e = 2 * norm_u_e / b_e.norm();
//		}
//
//		Ce(0, 0) += _configuration.sigma * h_e * norm_u_e;
//		Ce(1, 1) += _configuration.sigma * h_e * norm_u_e;
//
//		if (_propertiesConfiguration.gradient) {
//			matGradient.multiply(dND, temp, 1, 1);
//		}
//		if (_propertiesConfiguration.flux) {
//			matFlux.multiply(Ce, dND * temp, 1, 1);
//		}
//	}
//
//	if (_propertiesConfiguration.gradient) {
//		solution[offset + SolutionIndex::GRADIENT]->data[e->domains().front()].push_back(matGradient(0, 0) / e->gaussePoints());
//		solution[offset + SolutionIndex::GRADIENT]->data[e->domains().front()].push_back(matGradient(1, 0) / e->gaussePoints());
//	}
//
//	if (_propertiesConfiguration.flux) {
//		solution[offset + SolutionIndex::FLUX]->data[e->domains().front()].push_back(matFlux(0, 0) / e->gaussePoints());
//		solution[offset + SolutionIndex::FLUX]->data[e->domains().front()].push_back(matFlux(1, 0) / e->gaussePoints());
//	}
}

void HeatTransfer2D::processSolution(const Step &step)
{
	bool postProcess = false;
//	if (_propertiesConfiguration.gradient) {
//		postProcess = true;
//		if (_instance->solutions[offset + SolutionIndex::GRADIENT] == NULL) {
//			_instance->solutions[offset + SolutionIndex::GRADIENT] = new Solution(*_mesh, "gradient", ElementType::ELEMENTS, { Property::GRADIENT_X, Property::GRADIENT_Y });
//		}
//		for (size_t p = 0; p < _mesh->parts(); p++) {
//			_instance->solutions[offset + SolutionIndex::GRADIENT]->data[p].clear();
//		}
//	}
//	if (_propertiesConfiguration.flux) {
//		postProcess = true;
//		if (_instance->solutions[offset + SolutionIndex::FLUX] == NULL) {
//			_instance->solutions[offset + SolutionIndex::FLUX] = new Solution(*_mesh, "flux", ElementType::ELEMENTS, { Property::FLUX_X, Property::FLUX_Y });
//		}
//		for (size_t p = 0; p < _mesh->parts(); p++) {
//			_instance->solutions[offset + SolutionIndex::FLUX]->data[p].clear();
//		}
//	}
//
//	bool phase_change = false;
//	for (size_t m = 0; m < _mesh->materials().size(); m++) {
//		phase_change |= _mesh->materials()[m]->phase_change;
//	}
//
//	if (_propertiesConfiguration.phase && phase_change) {
//		postProcess = true;
//		if (_instance->solutions[offset + SolutionIndex::PHASE] == NULL) {
//			_instance->solutions[offset + SolutionIndex::PHASE] = new Solution(*_mesh, "phase", ElementType::NODES, { Property::PHASE });
//		}
//		for (size_t p = 0; p < _mesh->parts(); p++) {
//			_instance->solutions[offset + SolutionIndex::PHASE]->data[p].clear();
//		}
//	}
//
//	if (_propertiesConfiguration.latent_heat && phase_change) {
//		postProcess = true;
//		if (_instance->solutions[offset + SolutionIndex::LATENT_HEAT] == NULL) {
//			_instance->solutions[offset + SolutionIndex::LATENT_HEAT] = new Solution(*_mesh, "latent_heat", ElementType::NODES, { Property::LATENT_HEAT });
//		}
//		for (size_t p = 0; p < _mesh->parts(); p++) {
//			_instance->solutions[offset + SolutionIndex::LATENT_HEAT]->data[p].clear();
//		}
//	}
//
//	if (postProcess) {
//		#pragma omp parallel for
//		for (size_t p = 0; p < _mesh->parts(); p++) {
//			for (eslocal e = _mesh->getPartition()[p]; e < _mesh->getPartition()[p + 1]; e++) {
//				postProcessElement(step, _mesh->elements()[e], _instance->solutions);
//			}
//		}
//	}
}

