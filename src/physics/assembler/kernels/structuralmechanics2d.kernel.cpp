
#include "physics/assembler/dataholder.h"
#include "esinfo/time.h"
#include "esinfo/ecfinfo.h"
#include "structuralmechanics2d.kernel.h"

#include "basis/containers/point.h"
#include "basis/matrices/denseMatrix.h"
#include "basis/evaluator/evaluator.h"
#include "basis/logging/logging.h"
#include "config/ecf/root.h"
#include "mesh/elements/element.h"

using namespace espreso;

using namespace espreso;

void StructuralMechanics2DKernel::assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, DenseMatrix &K) const
{
	double Ex = 0, Ey = 0, mi = 0;
	Point p(coordinates[0], coordinates[1], 0);

	switch (mat->linear_elastic_properties.model) {

	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
		Ex = Ey = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate(p, time, temp);
		mi = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate(p, time, temp);
		break;

	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
		ESINFO(ERROR) << "Implement ANISOTROPIC MATERIAL";
		break;

	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
		Ex = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate(p, time, temp);
		Ey = mat->linear_elastic_properties.young_modulus.get(1, 1).evaluator->evaluate(p, time, temp);
		mi = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate(p, time, temp);
		break;

	default:
		ESINFO(ERROR) << "Linear elasticity 2D not supports set material model";
	}

	switch (mat->linear_elastic_properties.model) {

	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
	{

		switch (info::ecf->structural_mechanics_2d.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
		{
			double k = Ex * (1 - mi) / ((1 + mi) * (1 - 2 * mi));
			K(node, 0) = k * 1;
			K(node, 1) = k * 1;
			K(node, 2) = k * ((1 - 2 * mi) / (2 * (1 - mi)));
			K(node, 3) = k * (mi / (1 - mi));
			K(node, 4) = 0;
			K(node, 5) = 0;
			K(node, 6) = k * (mi / (1 - mi));
			K(node, 7) = 0;
			K(node, 8) = 0;
			return;
		}

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
		{
			double k = Ex / (1 - mi * mi);
			K(node, 0) = k * 1;
			K(node, 1) = k * 1;
			K(node, 2) = k * ((1 -  mi) / 2);
			K(node, 3) = k * mi;
			K(node, 4) = 0;
			K(node, 5) = 0;
			K(node, 6) = k * mi;
			K(node, 7) = 0;
			K(node, 8) = 0;
			return;
		}

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
		{
			K.resize(K.rows(), 16);
			double k = Ex * (1 - mi) / ((1 + mi) * (1 - 2 * mi));
			K(node,  0) = k * 1;
			K(node,  1) = k * 1;
			K(node,  2) = k * ((1 - 2 * mi) / (2 * (1 - mi)));
			K(node,  3) = k * 1;

			K(node,  4) = k * (mi / (1 - mi));
			K(node,  5) = 0;
			K(node,  6) = k * (mi / (1 - mi));
			K(node,  7) = 0;
			K(node,  8) = k * (mi / (1 - mi));
			K(node,  9) = 0;

			K(node, 10) = k * (mi / (1 - mi));
			K(node, 11) = 0;
			K(node, 12) = 0;
			K(node, 13) = k * (mi / (1 - mi));
			K(node, 14) = k * (mi / (1 - mi));
			K(node, 15) = 0;
			return;
		}
		}
		break;
	}

	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
	{
		ESINFO(ERROR) << "IMPLEMENT: MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC";
		return;
	}

	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
	{
		ESINFO(ERROR) << "IMPLEMENT: MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC";
		return;
	}

	default:
		ESINFO(ERROR) << "Structural mechanics 2D not supports set material model";
	}
}

void StructuralMechanics2DKernel::processElement(Matrices matrices, const SolverParameters &parameters, const ElementIterator &iterator, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const
{
	esint size = iterator.element->nodes;

	const std::vector<DenseMatrix> &N = *(iterator.element->N);
	const std::vector<DenseMatrix> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	DenseMatrix Ce(4, 4), XY(1, 2), coordinates(size, 2), J, invJ(2, 2), dND, B, precision, rhsT;
	DenseMatrix K(size, 9), TE(size, 2), thickness(size, 1), inertia(size, 2), dens(size, 1);
	DenseMatrix gpK(size, 9), gpTE(1, 2), gpThickness(1, 1), gpInertia(1, 2), gpDens(1, 1);
	double detJ, CP = 1, te;
	Point center;

	for (esint n = 0; n < size; n++) {
		inertia(n, 0) = iterator.acceleration[2 * n + 0];
		inertia(n, 1) = iterator.acceleration[2 * n + 1];
		coordinates(n, 0) = iterator.coordinates[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates[2 * n + 1];
		iterator.material->density.evaluator->evalVector(1, 2, iterator.coordinates + 2 * n, iterator.temperature + n, time::current, &dens(n, 0));

		thickness(n, 0) = iterator.thickness[n];
		switch (iterator.material->linear_elastic_properties.model) {
		case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
			iterator.material->linear_elastic_properties.thermal_expansion.get(0, 0).evaluator->evalVector(1, 2, iterator.coordinates, iterator.temperature, time::current, &te);
			TE(n, 0) = TE(n, 1) = (iterator.temperature[n] - iterator.initialTemperature[n]) * te;
			break;
		case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
			iterator.material->linear_elastic_properties.thermal_expansion.get(0, 0).evaluator->evalVector(1, 2, iterator.coordinates, iterator.temperature, time::current, &te);
			TE(n, 0) = (iterator.temperature[n] - iterator.initialTemperature[n]) * te;
			iterator.material->linear_elastic_properties.thermal_expansion.get(1, 1).evaluator->evalVector(1, 2, iterator.coordinates, iterator.temperature, time::current, &te);
			TE(n, 1) = (iterator.temperature[n] - iterator.initialTemperature[n]) * te;
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Invalid LINEAR ELASTIC model.";
		}
		assembleMaterialMatrix(n, iterator.coordinates, iterator.material, time::current, iterator.temperature[n], K);
	}

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if (matrices & (Matrices::K | Matrices::R)) {
		Ke.resize(2 * size, 2 * size);
		Ke = 0;
	}
	if (matrices & Matrices::M) {
		Me.resize(2 * size, 2 * size);
		Me = 0;
	}
	if (matrices & Matrices::R) {
		Re.resize(2 * size, 1);
		Re = 0;
	}
	if (matrices & Matrices::f) {
		fe.resize(2 * size, 1);
		fe = 0;
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J.values());
		inverse2x2(J.values(), invJ.values(), detJ);

		gpThickness.multiply(N[gp], thickness);
		gpK.multiply(N[gp], K);
		dND.multiply(invJ, dN[gp]);
		gpDens.multiply(N[gp], dens);

		if (matrices & Matrices::f) {
			gpTE.multiply(N[gp], TE);
			gpInertia.multiply(N[gp], inertia);
			XY.multiply(N[gp], coordinates);
		}

		if (matrices & Matrices::M) {
			Me.multiply(N[gp], N[gp], gpDens(0, 0) * detJ * weighFactor[gp] * CP, 1, true);
		}

		switch (info::ecf->structural_mechanics_2d.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:

			Ce.resize(3, 3);
			Ce(0, 0) = gpK(0, 0);
			Ce(1, 1) = gpK(0, 1);
			Ce(2, 2) = gpK(0, 2);
			Ce(0, 1) = gpK(0, 3);
			Ce(0, 2) = gpK(0, 4);
			Ce(1, 2) = gpK(0, 5);
			Ce(1, 0) = gpK(0, 6);
			Ce(2, 0) = gpK(0, 7);
			Ce(2, 1) = gpK(0, 8);

			B.resize(Ce.rows(), 2 * size);
			distribute3x2(B.values(), dND.values(), dND.rows(), dND.columns());

			if (matrices & (Matrices::K | Matrices::R)) {
				Ke.multiply(B, Ce * B, detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
			}

			if (matrices & Matrices::f) {
				precision.resize(Ce.rows(), 1);
				precision(0, 0) = gpTE(0, 1);
				precision(1, 0) = gpTE(0, 1);
				precision(2, 0) = 0;
				rhsT.multiply(B, Ce * precision, detJ * weighFactor[gp] * gpThickness(0, 0), 0, true, false);
				for (esint i = 0; i < 2 * size; i++) {
					fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % size) * gpInertia(0, i / size);
					fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * gpThickness(0, 0) * N[gp](0, i % size) * XY(0, i / size) * pow(iterator.angularVelocity[2], 2);
					fe(i, 0) += rhsT(i, 0);
				}
			}
			break;

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:

			Ce.resize(4, 4);
			Ce(0, 0) = gpK(0,  0);
			Ce(1, 1) = gpK(0,  1);
			Ce(2, 2) = gpK(0,  2);
			Ce(3, 3) = gpK(0,  3);
			Ce(0, 1) = gpK(0,  4);
			Ce(0, 2) = gpK(0,  5);
			Ce(0, 3) = gpK(0,  6);
			Ce(1, 2) = gpK(0,  7);
			Ce(1, 3) = gpK(0,  8);
			Ce(2, 3) = gpK(0,  9);
			Ce(1, 0) = gpK(0, 10);
			Ce(2, 0) = gpK(0, 11);
			Ce(2, 1) = gpK(0, 12);
			Ce(3, 0) = gpK(0, 13);
			Ce(3, 1) = gpK(0, 14);
			Ce(3, 2) = gpK(0, 15);

			B.resize(Ce.rows(), 2 * size);
			distribute4x2(B.values(), dND.values(), dND.rows(), dND.columns());
			for(size_t i = 0; i < N[gp].columns(); i++) {
				B(2, i) = N[gp](0, i) / XY(0, 0);
			}

			if (matrices & (Matrices::K | Matrices::R)) {
				Ke.multiply(B, Ce * B, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 1, true);
			}

			if (matrices & Matrices::f) {
				precision.resize(Ce.rows(), 1);
				precision(0, 0) = gpTE(0, 0);
				precision(1, 0) = gpTE(0, 1);
				precision(2, 0) = precision(3, 0) = 0;
				rhsT.multiply(B, Ce * precision, detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0), 0, true, false);
				for (esint i = 0; i < 2 * size; i++) {
					fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % size) * gpInertia(0, i / size);
					fe(i, 0) += rhsT(i, 0);
				}
				for (esint i = 0; i < size; i++) {
					fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % size) * XY(0, 0) * pow(iterator.angularVelocity[1], 2);
					fe(size + i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % size) * XY(0, 1) * pow(iterator.angularVelocity[1], 2);
				}
			}
			break;
		}
	}
}

void StructuralMechanics2DKernel::processEdge(Matrices matrices, const SolverParameters &parameters, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const
{
	Ke.resize(0, 0);
	if (iterator.normalPressure == NULL) {
		fe.resize(0, 0);
		return;
	}
	if (!(matrices & (Matrices::K | Matrices::f))) {
		fe.resize(0, 0);
		return;
	}

	esint size = iterator.element->nodes;

	const std::vector<DenseMatrix> &N = *(iterator.element->N);
	const std::vector<DenseMatrix> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	DenseMatrix coordinates(size, 2), dND(1, 2), P(size, 1), normal(2, 1), matThickness(size, 1), XY(1, 2);
	DenseMatrix gpP(1, 1), gpQ(1, 2), gpThickness(1, 1);

	fe.resize(0, 0);

	if (matrices & Matrices::f) {
		fe.resize(2 * size, 1);
		fe = 0;
	}

	for (esint n = 0; n < size; n++) {
		coordinates(n, 0) = iterator.coordinates[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates[2 * n + 1];
		P(n, 0) = iterator.normalPressure[n];
		matThickness(n, 0) = iterator.thickness[n];
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		Point n(-dND(0, 1), dND(0, 0), 0);
		// e->rotateOutside(e->parentElements()[0], _mesh->coordinates(), n);
		normal(0, 0) = n.x / J;
		normal(1, 0) = n.y / J;
		gpP.multiply(N[gp], P);
		gpQ.multiply(normal, gpP);
		gpThickness.multiply(N[gp], matThickness);

		switch (info::ecf->structural_mechanics_2d.element_behaviour) {

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			gpThickness(0, 0) = 1;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			for (esint i = 0; i < 2 * size; i++) {
				fe(i, 0) += gpThickness(0, 0) * J * weighFactor[gp] * N[gp](0, i % size) * gpQ(0, i / size);
			}
			break;

		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			XY.multiply(N[gp], coordinates);
			for (esint i = 0; i < 2 * size; i++) {
				fe(i, 0) += gpThickness(0, 0) * J * weighFactor[gp] * 2 * M_PI * XY(0, 0) * N[gp](0, i % size) * gpQ(0, i / size);
			}
			break;
		}
	}
}

void StructuralMechanics2DKernel::processSolution(const SolutionIterator &iterator)
{

}





