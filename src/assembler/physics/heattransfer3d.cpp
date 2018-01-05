
#include "heattransfer3d.h"

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

#ifdef BEM4I
#include "esbem.h"
#endif

// TODO: create file fith constants
#define CONST_Stefan_Boltzmann 5.6703e-8

using namespace espreso;

HeatTransfer3D::HeatTransfer3D(Mesh *mesh, Instance *instance, Step *step, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration)
: Physics("HEAT TRANSFER 3D", mesh, instance, step, &configuration), HeatTransfer(configuration, propertiesConfiguration)
{
	if (_propertiesConfiguration.gradient) {
		_gradient = _mesh->elements->appendData({ "GRADIENT", "GRADIENT_X", "GRADIENT_Y", "GRADIENT_Z" });
	}

	if (_propertiesConfiguration.flux) {
		_flux = _mesh->elements->appendData({ "FLUX", "FLUX_X", "FLUX_Y", "FLUX_Z" });
	}
}

void HeatTransfer3D::processBEM(eslocal domain, Matrices matrices)
{
	_instance->K[domain].rows = _mesh->nodes->dintervals[domain].back().DOFOffset;
	_instance->K[domain].cols = _instance->K[domain].rows;
	_instance->K[domain].nnz  = _instance->K[domain].rows * _instance->K[domain].cols;
	_instance->K[domain].type = 'G';
	_instance->K[domain].dense_values.resize(_instance->K[domain].nnz);
	_instance->K[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;

#ifndef BEM4I
	ESINFO(GLOBAL_ERROR) << "BEM4I is not linked!";
#else
	std::vector<eslocal> elements;
	std::vector<double> coordinates;

	bem4i::getLaplaceSteklovPoincare(
			_instance->K[domain].dense_values.data(),
			_mesh->nodes->dintervals[domain].back().DOFOffset,
			coordinates.data(),
			(eslocal)(elements.size() / 3),
			elements.data(),
			0,
			3, 3,
			&_BEMData[domain],
			0);
#endif

}

void HeatTransfer3D::assembleMaterialMatrix(eslocal eindex, eslocal node, const Point &p, const MaterialBaseConfiguration *mat, double phase, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const
{
	auto d2r = [] (double degree) -> double {
		return M_PI * degree / 180;
	};

	Point sin, cos;
	switch (mat->coordinate_system.type) {
	case CoordinateSystemConfiguration::TYPE::CARTESIAN:

		cos.x = std::cos(d2r(mat->coordinate_system.rotation_x.evaluator->evaluate(p, _step->currentTime, temp)));
		cos.y = std::cos(d2r(mat->coordinate_system.rotation_y.evaluator->evaluate(p, _step->currentTime, temp)));
		cos.z = std::cos(d2r(mat->coordinate_system.rotation_z.evaluator->evaluate(p, _step->currentTime, temp)));

		sin.x = std::sin(d2r(mat->coordinate_system.rotation_x.evaluator->evaluate(p, _step->currentTime, temp)));
		sin.y = std::sin(d2r(mat->coordinate_system.rotation_y.evaluator->evaluate(p, _step->currentTime, temp)));
		sin.z = std::sin(d2r(mat->coordinate_system.rotation_z.evaluator->evaluate(p, _step->currentTime, temp)));

		break;

	case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: {

		Point origin(
				mat->coordinate_system.center_x.evaluator->evaluate(p, _step->currentTime, temp),
				mat->coordinate_system.center_y.evaluator->evaluate(p, _step->currentTime, temp),
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
				mat->coordinate_system.center_x.evaluator->evaluate(p, _step->currentTime, temp),
				mat->coordinate_system.center_y.evaluator->evaluate(p, _step->currentTime, temp),
				mat->coordinate_system.center_z.evaluator->evaluate(p, _step->currentTime, temp));

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
				expression.evaluator->evaluate(p, _step->currentTime, temp + h) -
				expression.evaluator->evaluate(p, _step->currentTime, temp - h)
				) / (2 * h);
	};

	switch (mat->thermal_conductivity.model) {
	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		C(0, 0) = C(1, 1) = C(2, 2) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = _CD(1, 1) = _CD(2, 2) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(0, 1) = _CD(0, 2) = _CD(1, 0) = _CD(1, 2) = _CD(2, 0) = _CD(2, 1) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		C(2, 2) = mat->thermal_conductivity.values.get(2, 2).evaluator->evaluate(p, _step->currentTime, temp);
		C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(2, 2) = derivation(mat->thermal_conductivity.values.get(2, 2), temp / 1e4);
			_CD(0, 1) = _CD(0, 2) = _CD(1, 0) = _CD(1, 2) = _CD(2, 0) = _CD(2, 1) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		C(2, 2) = mat->thermal_conductivity.values.get(2, 2).evaluator->evaluate(p, _step->currentTime, temp);
		C(0, 1) = C(1, 0) = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate(p, _step->currentTime, temp);
		C(0, 2) = C(2, 0) = mat->thermal_conductivity.values.get(0, 2).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 2) = C(2, 1) = mat->thermal_conductivity.values.get(1, 2).evaluator->evaluate(p, _step->currentTime, temp);
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
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate(p, _step->currentTime, temp);
		C(2, 2) = mat->thermal_conductivity.values.get(2, 2).evaluator->evaluate(p, _step->currentTime, temp);
		C(0, 1) = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate(p, _step->currentTime, temp);
		C(0, 2) = mat->thermal_conductivity.values.get(0, 2).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 0) = mat->thermal_conductivity.values.get(1, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(1, 2) = mat->thermal_conductivity.values.get(1, 2).evaluator->evaluate(p, _step->currentTime, temp);
		C(2, 0) = mat->thermal_conductivity.values.get(2, 0).evaluator->evaluate(p, _step->currentTime, temp);
		C(2, 1) = mat->thermal_conductivity.values.get(2, 1).evaluator->evaluate(p, _step->currentTime, temp);
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

void HeatTransfer3D::processElement(eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	auto nodes = _mesh->elements->nodes->cbegin() + eindex;
	auto epointer = _mesh->elements->epointers->datatarray()[eindex];
	const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[domain];
	const ECFExpressionVector *translation_motion = NULL;
	Evaluator *heat_source = NULL;
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

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	bool CAU = _configuration.stabilization == HeatTransferConfiguration::STABILIZATION::CAU;
	bool tangentCorrection = (matrices & Matrices::K) && _step->tangentMatrixCorrection;

	DenseMatrix Ce(3, 3), coordinates(nodes->size(), 3), J(3, 3), invJ(3, 3), dND;
	double detJ, temp;
	DenseMatrix f(nodes->size(), 1);
	DenseMatrix U(nodes->size(), 3);
	DenseMatrix m(nodes->size(), 1);
	DenseMatrix T(nodes->size(), 1);
	DenseMatrix K(nodes->size(), 9);
	DenseMatrix gpK(1, 9), gpM(1, 1);
	DenseMatrix tangentK, BT, BTN, gpCD, CD, CDBTN, CDe;

	const MaterialConfiguration* material = _mesh->materials[_mesh->elements->material->datatarray()[eindex]];

	const MaterialBaseConfiguration *phase1, *phase2;
	if (material->phase_change) {
		phase1 = &material->phases.find(1)->second;
		phase2 = &material->phases.find(2)->second;
	}

	if (tangentCorrection) {
		CD.resize(nodes->size(), 9);
		CDe.resize(3, 3);
	}

	for (size_t n = 0; n < nodes->size(); n++) {
		auto it = std::lower_bound(intervals.begin(), intervals.end(), nodes->at(n), [] (const DomainInterval &interval, eslocal node) { return interval.end < node; });
		temp = (*_temperature->decomposedData)[domain][it->DOFOffset + nodes->at(n) - it->begin];
		const Point &p = _mesh->nodes->coordinates->datatarray()[nodes->at(n)];
		T(n, 0) = temp;
		coordinates(n, 0) = p.x;
		coordinates(n, 1) = p.y;
		coordinates(n, 2) = p.z;
		if (material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, material->phase_change_temperature - material->transition_interval / 2, material->phase_change_temperature + material->transition_interval / 2, temp, material->smooth_step_order);
			assembleMaterialMatrix(eindex, n, p, phase1, phase, temp, K, CD, tangentCorrection);
			assembleMaterialMatrix(eindex, n, p, phase2, (1 - phase), temp, K, CD, tangentCorrection);
			m(n, 0) =
					(    phase  * phase1->density.evaluator->evaluate(p, _step->currentTime, temp) +
					(1 - phase) * phase2->density.evaluator->evaluate(p, _step->currentTime, temp)) *

					(    phase  * phase1->heat_capacity.evaluator->evaluate(p, _step->currentTime, temp) +
					(1 - phase) * phase2->heat_capacity.evaluator->evaluate(p, _step->currentTime, temp) +
					material->latent_heat * derivation);
		} else {
			assembleMaterialMatrix(eindex, n, p, material, 1, temp, K, CD, tangentCorrection);
			m(n, 0) =
					material->density.evaluator->evaluate(p, _step->currentTime, temp) *
					material->heat_capacity.evaluator->evaluate(p, _step->currentTime, temp);
		}

		if (translation_motion) {
			U(n, 0) = translation_motion->x.evaluator->evaluate(p, _step->currentTime, temp) * m(n, 0);
			U(n, 1) = translation_motion->y.evaluator->evaluate(p, _step->currentTime, temp) * m(n, 0);
			U(n, 2) = translation_motion->z.evaluator->evaluate(p, _step->currentTime, temp) * m(n, 0);
		}
		if (heat_source) {
			f(n, 0) = heat_source->evaluate(p, _step->currentTime, temp);
		}
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

	DenseMatrix u(1, 3), v(1, 3), re(1, nodes->size());
	double normGradN = 0;

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

		DenseMatrix b_e(1, nodes->size()), b_e_c(1, nodes->size());
		b_e.multiply(u, dND, 1, 0);

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
		double h_e = 0, tau_e = 0, konst = 0;
		double C_e = 0;

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

		Ce(0, 0) += _configuration.sigma * h_e * norm_u_e;
		Ce(1, 1) += _configuration.sigma * h_e * norm_u_e;
		Ce(2, 2) += _configuration.sigma * h_e * norm_u_e;

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
			Ke.multiply(dND, Ce * dND, detJ * weighFactor[gp], 1, true);
			Ke.multiply(N[gp], b_e, detJ * weighFactor[gp], 1, true);
			if (konst * weighFactor[gp] * detJ != 0) {
				Ke.multiply(b_e, b_e, konst * weighFactor[gp] * detJ, 1, true);
			}
			if (CAU) {
				Ke.multiply(dND, dND, C_e * weighFactor[gp] * detJ, 1, true);
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
}

void HeatTransfer3D::processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	const ConvectionConfiguration *convection = NULL;
	const RadiationConfiguration *radiation = NULL;
	const Evaluator *heatFlow = NULL, *heatFlux = NULL;

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

	auto nodes = region->elements->cbegin() + findex;
	auto epointer = region->epointers->datatarray()[findex];
	const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[domain];

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	DenseMatrix coordinates(nodes->size(), 3), dND(1, 3), q(nodes->size(), 1), htc(nodes->size(), 1), flow(nodes->size(), 1), emiss(nodes->size(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpFlow(1, 1), gpEmiss(1, 1);

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
		coordinates(n, 2) = p.z;
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
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		double J = va.norm();

		gpQ.multiply(N[gp], q);
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

void HeatTransfer3D::processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
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

void HeatTransfer3D::processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
}

void HeatTransfer3D::postProcessElement(eslocal domain, eslocal eindex)
{
	auto nodes = _mesh->elements->nodes->cbegin() + eindex;
	auto epointer = _mesh->elements->epointers->datatarray()[eindex];
	const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[domain];
	const ECFExpressionVector *translation_motion = NULL;
	for (auto it = _configuration.load_steps_settings.at(_step->step + 1).translation_motions.begin(); it != _configuration.load_steps_settings.at(_step->step + 1).translation_motions.end(); ++it) {
		ElementsRegionStore *region = _mesh->eregion(it->first);
		if (std::binary_search(region->elements->datatarray().cbegin(), region->elements->datatarray().cend(), eindex)) {
			translation_motion = &it->second;
			break;
		}
	}

	const std::vector<DenseMatrix> &N = *(epointer->N);
	const std::vector<DenseMatrix> &dN = *(epointer->dN);
	const std::vector<double> &weighFactor = *(epointer->weighFactor);

	DenseMatrix Ce(3, 3), coordinates, J(3, 3), invJ(3, 3), dND, temp(nodes->size(), 1);
	double detJ, m, norm_u_e, h_e;
	DenseMatrix U(nodes->size(), 3), K(nodes->size(), 9), gpK(1, 9), CD;
	DenseMatrix u(1, 3), matFlux(3, 1), matGradient(3, 1);

	const MaterialConfiguration* material = _mesh->materials[_mesh->elements->material->datatarray()[eindex]];

	const MaterialBaseConfiguration *phase1, *phase2;
	if (material->phase_change) {
		phase1 = &material->phases.find(1)->second;
		phase2 = &material->phases.find(2)->second;
	}

	coordinates.resize(nodes->size(), 3);

	for (size_t i = 0; i < nodes->size(); i++) {
		auto it = std::lower_bound(intervals.begin(), intervals.end(), nodes->at(i), [] (const DomainInterval &interval, eslocal node) { return interval.end < node; });
		temp(i, 0) = (*_temperature->decomposedData)[domain][it->DOFOffset + nodes->at(i) - it->begin];
		const Point &p = _mesh->nodes->coordinates->datatarray()[nodes->at(i)];
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		coordinates(i, 2) = p.z;
		if (material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, material->phase_change_temperature - material->transition_interval / 2, material->phase_change_temperature + material->transition_interval / 2, temp(i, 0), material->smooth_step_order);
			assembleMaterialMatrix(eindex, i, p, phase1, phase, temp(i, 0), K, CD, false);
			assembleMaterialMatrix(eindex, i, p, phase2, (1 - phase), temp(i, 0), K, CD, false);
			m =
					(    phase  * phase1->density.evaluator->evaluate(p, _step->currentTime, temp(i, 0)) +
					(1 - phase) * phase2->density.evaluator->evaluate(p, _step->currentTime, temp(i, 0))) *

					(    phase  * phase1->heat_capacity.evaluator->evaluate(p, _step->currentTime, temp(i, 0)) +
					(1 - phase) * phase2->heat_capacity.evaluator->evaluate(p, _step->currentTime, temp(i, 0)) +
					material->latent_heat * derivation);

			if (_phaseChange) {
				(*_phaseChange->decomposedData)[domain][it->DOFOffset + nodes->at(i) - it->begin] = phase;
				(*_latentHeat->decomposedData)[domain][it->DOFOffset + nodes->at(i) - it->begin] = material->latent_heat * derivation;
			}
		} else {
			assembleMaterialMatrix(eindex, i, p, material, 1, temp(i, 0), K, CD, false);
			m =
					material->density.evaluator->evaluate(p, _step->currentTime, temp(i, 0)) *
					material->heat_capacity.evaluator->evaluate(p, _step->currentTime, temp(i, 0));
		}

		if (translation_motion) {
			U(i, 0) = translation_motion->x.evaluator->evaluate(p, _step->currentTime, temp(i, 0));
			U(i, 1) = translation_motion->y.evaluator->evaluate(p, _step->currentTime, temp(i, 0));
			U(i, 2) = translation_motion->z.evaluator->evaluate(p, _step->currentTime, temp(i, 0));
		}
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
			DenseMatrix b_e(1, nodes->size());
			b_e.multiply(u, dND, 1, 0);
			h_e = 2 * norm_u_e / b_e.norm();
		}

		Ce(0, 0) += _configuration.sigma * h_e * norm_u_e;
		Ce(1, 1) += _configuration.sigma * h_e * norm_u_e;
		Ce(2, 2) += _configuration.sigma * h_e * norm_u_e;

		matGradient.multiply(dND, temp, 1, 1);
		matFlux.multiply(Ce, dND * temp, 1, 1);
	}

	if (_propertiesConfiguration.gradient) {
		(*_gradient->data)[3 * eindex + 0] = matGradient(0, 0) / N.size();
		(*_gradient->data)[3 * eindex + 1] = matGradient(1, 0) / N.size();
		(*_gradient->data)[3 * eindex + 2] = matGradient(2, 0) / N.size();
	}

	if (_propertiesConfiguration.flux) {
		(*_flux->data)[3 * eindex + 0] = matFlux(0, 0) / N.size();
		(*_flux->data)[3 * eindex + 1] = matFlux(1, 0) / N.size();
		(*_flux->data)[3 * eindex + 2] = matFlux(2, 0) / N.size();
	}
}

void HeatTransfer3D::processSolution()
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
