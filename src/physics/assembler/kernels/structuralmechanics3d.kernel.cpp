
#include "esinfo/timeinfo.h"
#include "esinfo/eslog.hpp"
#include "physics/assembler/dataholder.h"
#include "structuralmechanics3d.kernel.h"

#include "basis/containers/point.h"
#include "basis/matrices/denseMatrix.h"
#include "basis/evaluator/evaluator.h"
#include "config/ecf/physics/heattransfer.h"
#include "mesh/elements/element.h"

#include "basis/utilities/print.h"

using namespace espreso;

using namespace espreso;

void StructuralMechanics3DKernel::assembleLinearElasticMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, DenseMatrix &K) const
{
	double Ex, Ey, Ez, miXY, miXZ, miYZ, Gx, Gy, Gz;
	Point p(coordinates[0], coordinates[1], coordinates[2]);

	switch (mat->linear_elastic_properties.model) {

	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC: {
		Ex = Ey = Ez = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate(p, time, temp);
		miXY = miXZ = miYZ = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate(p, time, temp);

		double EE = Ex / ((1 + miXY) * (1 - 2 * miXY));

		K(node,  0) = EE * (1.0 - miXY);
		K(node,  1) = EE * (1.0 - miXY);
		K(node,  2) = EE * (1.0 - miXY);
		K(node,  3) = EE * (0.5 - miXY);
		K(node,  4) = EE * (0.5 - miXY);
		K(node,  5) = EE * (0.5 - miXY);

		K(node,  6) = EE * miXY;
		K(node,  7) = EE * miXY;
		K(node,  8) = 0;
		K(node,  9) = 0;
		K(node, 10) = 0;
		K(node, 11) = EE * miXY;
		K(node, 12) = 0;
		K(node, 13) = 0;
		K(node, 14) = 0;
		K(node, 15) = 0;
		K(node, 16) = 0;
		K(node, 17) = 0;
		K(node, 18) = 0;
		K(node, 19) = 0;
		K(node, 20) = 0;
		K(node, 21) = EE * miXY;
		K(node, 22) = EE * miXY;
		K(node, 23) = EE * miXY;
		K(node, 24) = 0;
		K(node, 25) = 0;
		K(node, 26) = 0;
		K(node, 27) = 0;
		K(node, 28) = 0;
		K(node, 29) = 0;
		K(node, 30) = 0;
		K(node, 31) = 0;
		K(node, 32) = 0;
		K(node, 33) = 0;
		K(node, 34) = 0;
		K(node, 35) = 0;
	} break;

	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC: {
		K(node,  0) = mat->linear_elastic_properties.anisotropic.get(0, 0).evaluator->evaluate(p, time, temp);
		K(node,  1) = mat->linear_elastic_properties.anisotropic.get(1, 1).evaluator->evaluate(p, time, temp);
		K(node,  2) = mat->linear_elastic_properties.anisotropic.get(2, 2).evaluator->evaluate(p, time, temp);
		K(node,  3) = mat->linear_elastic_properties.anisotropic.get(3, 3).evaluator->evaluate(p, time, temp);
		K(node,  4) = mat->linear_elastic_properties.anisotropic.get(4, 4).evaluator->evaluate(p, time, temp);
		K(node,  5) = mat->linear_elastic_properties.anisotropic.get(5, 5).evaluator->evaluate(p, time, temp);

		K(node,  6) = mat->linear_elastic_properties.anisotropic.get(0, 1).evaluator->evaluate(p, time, temp);
		K(node,  7) = mat->linear_elastic_properties.anisotropic.get(0, 2).evaluator->evaluate(p, time, temp);
		K(node,  8) = mat->linear_elastic_properties.anisotropic.get(0, 3).evaluator->evaluate(p, time, temp);
		K(node,  9) = mat->linear_elastic_properties.anisotropic.get(0, 4).evaluator->evaluate(p, time, temp);
		K(node, 10) = mat->linear_elastic_properties.anisotropic.get(0, 5).evaluator->evaluate(p, time, temp);
		K(node, 11) = mat->linear_elastic_properties.anisotropic.get(1, 2).evaluator->evaluate(p, time, temp);
		K(node, 12) = mat->linear_elastic_properties.anisotropic.get(1, 3).evaluator->evaluate(p, time, temp);
		K(node, 13) = mat->linear_elastic_properties.anisotropic.get(1, 4).evaluator->evaluate(p, time, temp);
		K(node, 14) = mat->linear_elastic_properties.anisotropic.get(1, 5).evaluator->evaluate(p, time, temp);
		K(node, 15) = mat->linear_elastic_properties.anisotropic.get(2, 3).evaluator->evaluate(p, time, temp);
		K(node, 16) = mat->linear_elastic_properties.anisotropic.get(2, 4).evaluator->evaluate(p, time, temp);
		K(node, 17) = mat->linear_elastic_properties.anisotropic.get(2, 5).evaluator->evaluate(p, time, temp);
		K(node, 18) = mat->linear_elastic_properties.anisotropic.get(3, 4).evaluator->evaluate(p, time, temp);
		K(node, 19) = mat->linear_elastic_properties.anisotropic.get(3, 5).evaluator->evaluate(p, time, temp);
		K(node, 20) = mat->linear_elastic_properties.anisotropic.get(4, 5).evaluator->evaluate(p, time, temp);

		K(node, 21) = K(node,  6);
		K(node, 22) = K(node,  7);
		K(node, 23) = K(node, 11);
		K(node, 24) = K(node,  8);
		K(node, 25) = K(node, 12);
		K(node, 26) = K(node, 15);
		K(node, 27) = K(node,  9);
		K(node, 28) = K(node, 13);
		K(node, 29) = K(node, 16);
		K(node, 30) = K(node, 18);
		K(node, 31) = K(node, 10);
		K(node, 32) = K(node, 14);
		K(node, 33) = K(node, 17);
		K(node, 34) = K(node, 19);
		K(node, 35) = K(node, 20);
	} break;

	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC: {
		Ex = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate(p, time, temp);
		Ey = mat->linear_elastic_properties.young_modulus.get(1, 1).evaluator->evaluate(p, time, temp);
		Ez = mat->linear_elastic_properties.young_modulus.get(2, 2).evaluator->evaluate(p, time, temp);

		miXY = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate(p, time, temp);
		miXZ = mat->linear_elastic_properties.poisson_ratio.get(1, 1).evaluator->evaluate(p, time, temp);
		miYZ = mat->linear_elastic_properties.poisson_ratio.get(2, 2).evaluator->evaluate(p, time, temp);

		Gx = mat->linear_elastic_properties.shear_modulus.get(0, 0).evaluator->evaluate(p, time, temp);
		Gy = mat->linear_elastic_properties.shear_modulus.get(1, 1).evaluator->evaluate(p, time, temp);
		Gz = mat->linear_elastic_properties.shear_modulus.get(2, 2).evaluator->evaluate(p, time, temp);

		double miYX = miXY * Ey / Ex;
		double miZY = miYZ * Ez / Ey;
		double miZX = miXZ * Ex / Ez;

		double ksi = 1 - (miXY * miYX + miYZ * miZY + miXZ * miZX) - (miXY * miYZ * miZX + miYX * miZY * miXZ);

		double dxx = Ex * (1 - miYZ * miZY) / ksi;
		double dxy = Ey * (miXY + miXZ * miZY) / ksi;
		double dxz = Ez * (miXZ + miYZ * miXY)  /ksi;
		double dyy = Ey * (1 - miXZ * miZX) / ksi;
		double dyz = Ez * (miYZ + miYX * miXZ) / ksi;
		double dzz = Ez * (1 - miYX * miXY) / ksi;

		K(node,  0) = dxx;
		K(node,  1) = dyy;
		K(node,  2) = dzz;
		K(node,  3) = Gx;
		K(node,  4) = Gz;
		K(node,  5) = Gy;

		K(node,  6) = dxy;
		K(node,  7) = dxz;
		K(node,  8) = 0;
		K(node,  9) = 0;
		K(node, 10) = 0;
		K(node, 11) = dyz;
		K(node, 12) = 0;
		K(node, 13) = 0;
		K(node, 14) = 0;
		K(node, 15) = 0;
		K(node, 16) = 0;
		K(node, 17) = 0;
		K(node, 18) = 0;
		K(node, 19) = 0;
		K(node, 20) = 0;
		K(node, 21) = dxy;
		K(node, 22) = dxz;
		K(node, 23) = dyz;
		K(node, 24) = 0;
		K(node, 25) = 0;
		K(node, 26) = 0;
		K(node, 27) = 0;
		K(node, 28) = 0;
		K(node, 29) = 0;
		K(node, 30) = 0;
		K(node, 31) = 0;
		K(node, 32) = 0;
		K(node, 33) = 0;
		K(node, 34) = 0;
		K(node, 35) = 0;
	} break;

	default:
		eslog::error("Structural mechanics 3D not supports set material model.\n");
	}
}

void StructuralMechanics3DKernel::assembleHyperElasticMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, DenseMatrix &K) const
{

}

void StructuralMechanics3DKernel::processElement(Matrices matrices, const SolverParameters &parameters, const ElementIterator &iterator, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Ce, DenseMatrix &Re, DenseMatrix &fe) const
{
	int size = iterator.element->nodes;

	const std::vector<DenseMatrix> &N = *(iterator.element->N);
	const std::vector<DenseMatrix> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	DenseMatrix C(6, 6), initCoordinates(size, 3), coordinates(size, 3), F(3, 3), J, JC, invJ(3, 3), dND, B, precision, rhsT;
	DenseMatrix eHat(3, 3), eVec(6, 1), sVec, S(9, 9), BL(6, 3 * size), GL(9, 3 * size);
	DenseMatrix K(size, 36), TE(size, 3), inertia(size, 3), dens(size, 1);
	DenseMatrix gpK(size, 36), gpTE(1, 3), gpInertia(1, 3), gpDens(1, 1);
	double detJ, CP = 1, te;

	for (int n = 0; n < size; n++) {
		inertia(n, 0) = iterator.acceleration[3 * n + 0];
		inertia(n, 1) = iterator.acceleration[3 * n + 1];
		inertia(n, 2) = iterator.acceleration[3 * n + 2];
		initCoordinates(n, 0) = iterator.coordinates[3 * n + 0];
		initCoordinates(n, 1) = iterator.coordinates[3 * n + 1];
		initCoordinates(n, 2) = iterator.coordinates[3 * n + 2];
		coordinates(n, 0) = initCoordinates(n, 0) + iterator.displacement[3 * n + 0];
		coordinates(n, 1) = initCoordinates(n, 1) + iterator.displacement[3 * n + 1];
		coordinates(n, 2) = initCoordinates(n, 2) + iterator.displacement[3 * n + 2];
		iterator.material->density.evaluator->evalVector(1, 3, iterator.coordinates + 3 * n, iterator.temperature + n, time::current, &dens(n, 0));

		switch (iterator.material->thermal_expansion.model) {
		case ThermalExpansionConfiguration::MODEL::ISOTROPIC:
			iterator.material->thermal_expansion.thermal_expansion.get(0, 0).evaluator->evalVector(1, 3, iterator.coordinates, iterator.temperature, time::current, &te);
			TE(n, 0) = TE(n, 1) = TE(n, 2) = (iterator.temperature[n] - iterator.initialTemperature[n]) * te;
			break;
		case ThermalExpansionConfiguration::MODEL::ORTHOTROPIC:
			iterator.material->thermal_expansion.thermal_expansion.get(0, 0).evaluator->evalVector(1, 3, iterator.coordinates, iterator.temperature, time::current, &te);
			TE(n, 0) = (iterator.temperature[n] - iterator.initialTemperature[n]) * te;
			iterator.material->thermal_expansion.thermal_expansion.get(1, 1).evaluator->evalVector(1, 3, iterator.coordinates, iterator.temperature, time::current, &te);
			TE(n, 1) = (iterator.temperature[n] - iterator.initialTemperature[n]) * te;
			iterator.material->thermal_expansion.thermal_expansion.get(2, 2).evaluator->evalVector(1, 3, iterator.coordinates, iterator.temperature, time::current, &te);
			TE(n, 2) = (iterator.temperature[n] - iterator.initialTemperature[n]) * te;
			break;
		}

		switch (iterator.material->material_model) {
		case MaterialBaseConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC:
			assembleLinearElasticMaterialMatrix(n, iterator.coordinates, iterator.material, time::current, iterator.temperature[n], K);
			break;
		case MaterialBaseConfiguration::MATERIAL_MODEL::HYPER_ELASTIC:
			assembleHyperElasticMaterialMatrix(n, iterator.coordinates, iterator.material, time::current, iterator.temperature[n], K);
			break;
		}
	}

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if (matrices & (Matrices::K | Matrices::R)) {
		Ke.resize(3 * size, 3 * size);
	}
	if (matrices & Matrices::M) {
		Me.resize(size, size);
		Me = 0;
	}
	if (matrices & Matrices::R) {
		Re.resize(3 * size, 1);
		Re = 0;
	}
	if (matrices & Matrices::f) {
		fe.resize(3 * size, 1);
		fe = 0;
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		J.multiply(dN[gp], initCoordinates);

		detJ = determinant3x3(J.values());
		inverse3x3(J.values(), invJ.values(), detJ);

		gpK.multiply(N[gp], K);
		dND.multiply(invJ, dN[gp]);
		gpDens.multiply(N[gp], dens);

		if (matrices & Matrices::f) {
			gpTE.multiply(N[gp], TE);
			gpInertia.multiply(N[gp], inertia);
		}

		if (matrices & Matrices::M) {
			Me.multiply(N[gp], N[gp], gpDens(0, 0) * detJ * weighFactor[gp] * CP, 1, true);
		}

		size_t k = 0;
		for (size_t i = 0; i < 6; i++) {
			C(i, i) = gpK(0, k++);
		}
		for (size_t i = 0; i < 6; i++) {
			for (size_t j = i + 1; j < 6; j++) {
				C(i, j) = gpK(0, k++);
			}
		}
		for (size_t i = 0; i < 6; i++) {
			for (size_t j = 0; j < i; j++) {
				C(i, j) = gpK(0, k++);
			}
		}

		if (iterator.largeDisplacement && (matrices & Matrices::R)) {
			JC.multiply(dN[gp], coordinates);
			F.multiply(JC, invJ, 1, 0, true);

			if (true) { // case LINEAR ELASTIC MODEL
				eHat.multiply(F, F, .5, 0, true, false);
				eHat(0, 0) -= .5;
				eHat(1, 1) -= .5;
				eHat(2, 2) -= .5;
				eVec(0, 0) = eHat(0, 0);
				eVec(1, 0) = eHat(1, 1);
				eVec(2, 0) = eHat(2, 2);
				eVec(3, 0) = 2 * eHat(0, 1);
				eVec(4, 0) = 2 * eHat(1, 2);
				eVec(5, 0) = 2 * eHat(0, 2);
				sVec.multiply(C, eVec);
				for (int i = 0; i < 3; i++) {
					S(0 + i * 3, 0 + i * 3) = sVec(0, 0);
					S(1 + i * 3, 1 + i * 3) = sVec(1, 0);
					S(2 + i * 3, 2 + i * 3) = sVec(2, 0);
					S(0 + i * 3, 1 + i * 3) = sVec(3, 0);
					S(1 + i * 3, 0 + i * 3) = sVec(3, 0);
					S(0 + i * 3, 2 + i * 3) = sVec(5, 0);
					S(2 + i * 3, 0 + i * 3) = sVec(5, 0);
					S(1 + i * 3, 2 + i * 3) = sVec(4, 0);
					S(2 + i * 3, 1 + i * 3) = sVec(4, 0);
				}
				distribute9x3(GL.values(), dND.values(), dND.rows(), dND.columns());
				for (int i = 0; i < size; i++) {
					for (int j = 0; j < 3; j++) {
						BL(0, i + j * size) = F(j, 0) * dND(0, i);
						BL(1, i + j * size) = F(j, 1) * dND(1, i);
						BL(2, i + j * size) = F(j, 2) * dND(2, i);
						BL(3, i + j * size) = F(j, 0) * dND(1, i) + F(j, 1) * dND(0, i);
						BL(4, i + j * size) = F(j, 1) * dND(2, i) + F(j, 2) * dND(1, i);
						BL(5, i + j * size) = F(j, 0) * dND(2, i) + F(j, 2) * dND(0, i);
					}
				}
			}
			if (matrices & Matrices::K) {
				Ke.multiply(BL, C * BL, detJ * weighFactor[gp], 1, true);
				Ke.multiply(GL, S * GL, detJ * weighFactor[gp], 1, true);
			}
			Re.multiply(BL, sVec, detJ * weighFactor[gp], 1, true);
		} else {
			B.resize(C.rows(), 3 * size);
			distribute6x3(B.values(), dND.values(), dND.rows(), dND.columns());
			if (matrices & Matrices::K) {
				Ke.multiply(B, C * B, detJ * weighFactor[gp], 1, true);
			}

			if (matrices & Matrices::f) {
				precision.resize(C.rows(), 1);
				precision(0, 0) = gpTE(0, 0);
				precision(1, 0) = gpTE(0, 1);
				precision(2, 0) = gpTE(0, 2);
				precision(3, 0) = precision(4, 0) = precision(5, 0) = 0;

				rhsT.multiply(B, C * precision, detJ * weighFactor[gp], 0, true, false);
				for (esint i = 0; i < 3 * size; i++) {
					fe(i, 0) += rhsT(i, 0);
				}
			}
		}

		for (esint i = 0; i < 3 * size; i++) {
			fe(i, 0) += gpDens(0, 0) * detJ * weighFactor[gp] * N[gp](0, i % size) * gpInertia(0, i / size);
		}
	}
}

void StructuralMechanics3DKernel::processFace(Matrices matrices, const SolverParameters &parameters, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const
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

	DenseMatrix coordinates(size, 3), dND(1, 3), P(size, 1), normal(1, 3);
	DenseMatrix gpP(1, 1), gpQ(1, 3);

	esint Ksize = 3 * size;
	Ke.resize(0, 0);
	fe.resize(0, 0);

	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	for (esint n = 0; n < size; n++) {
		coordinates(n, 0) = iterator.coordinates[3 * n + 0];
		coordinates(n, 1) = iterator.coordinates[3 * n + 1];
		coordinates(n, 2) = iterator.coordinates[3 * n + 2];
		P(n, 0) += iterator.normalPressure[n];
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		// e->rotateOutside(e->parentElements()[0], _mesh->coordinates(), va);
		double J = va.norm();
		normal(0, 0) = va.x / va.norm();
		normal(0, 1) = va.y / va.norm();
		normal(0, 2) = va.z / va.norm();

		gpP.multiply(N[gp], P);
		gpQ.multiply(normal, gpP, 1, 0, true);

		for (esint i = 0; i < Ksize; i++) {
			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % size) * gpQ(0, i / size);
		}
	}
}

void StructuralMechanics3DKernel::processEdge(Matrices matrices, const SolverParameters &parameters, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const
{

}

void StructuralMechanics3DKernel::processSolution(const SolutionIterator &iterator)
{

}

