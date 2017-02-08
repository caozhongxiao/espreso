
#include "../../old_physics/elasticity3d/assembler.h"

#include "mkl.h"

#include "../../../basis/matrices/denseMatrix.h"
#include "../../../basis/matrices/sparseVVPMatrix.h"
#include "../../../basis/matrices/sparseCSRMatrix.h"
#include "../../../configuration/physics/linearelasticity3d.h"
#include "../../../solver/generic/SparseMatrix.h"
#include "../../../solver/specific/sparsesolvers.h"

#include "../../../mesh/elements/element.h"
#include "../../../mesh/settings/evaluator.h"
#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/structures/material.h"

#include "../../../mesh/elements/plane/square4.h"
#include "../../../mesh/elements/plane/square8.h"
#include "../../../mesh/elements/plane/triangle3.h"
#include "../../../mesh/elements/plane/triangle6.h"

#include "../../../mesh/elements/volume/hexahedron20.h"
#include "../../../mesh/elements/volume/hexahedron8.h"
#include "../../../mesh/elements/volume/prisma15.h"
#include "../../../mesh/elements/volume/prisma6.h"
#include "../../../mesh/elements/volume/pyramid13.h"
#include "../../../mesh/elements/volume/pyramid5.h"
#include "../../../mesh/elements/volume/tetrahedron10.h"
#include "../../../mesh/elements/volume/tetrahedron4.h"

#include "../../constraints/equalityconstraints.h"
#include "../../constraints/inequalityconstraints.h"

#include "../../../output/vtk/vtk.h"

namespace espreso {

std::vector<Property> Elasticity3D::elementDOFs;
std::vector<Property> Elasticity3D::faceDOFs;
std::vector<Property> Elasticity3D::edgeDOFs;
std::vector<Property> Elasticity3D::pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };
std::vector<Property> Elasticity3D::midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };

Elasticity3D::Elasticity3D(Mesh &mesh, Constraints &constraints, const LinearElasticity3DConfiguration &configuration)
: OldPhysics(
		mesh, constraints, configuration.espreso,
		MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
		elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs),
  _configuration(configuration) {};

void Elasticity3D::assembleStiffnessMatrices()
{
	ESINFO(PROGRESS2) << "Assemble matrices K, kernels, and RHS.";
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh.parts(); p++) {
		composeSubdomain(p);
		K[p].mtype = mtype;
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}

void Elasticity3D::prepareMeshStructures()
{
	Hexahedron8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Hexahedron20::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Tetrahedron4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Tetrahedron10::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Prisma6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Prisma15::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Pyramid5::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Pyramid13::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	Square4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Square8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle3::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	matrixSize = _mesh.assignUniformDOFsIndicesToNodes(matrixSize, pointDOFs);
	_mesh.computeNodesDOFsCounters(pointDOFs);

	if (_solverConfiguration.regularization == REGULARIZATION::FIX_POINTS) {
		_mesh.computeFixPoints(8);
	}

	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			_mesh.computeVolumeCorners(1, true, true, false);
			break;
		case B0_TYPE::KERNELS:
			_mesh.computeFacesSharedByDomains();
			break;
		case B0_TYPE::COMBINED:
			_mesh.computeFacesSharedByDomains();
			if (!_mesh.corners().size()) {
				_mesh.computeEdgesOnBordersOfFacesSharedByDomains();
				_mesh.computeCornersOnEdges(1, true, true);
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented B0";
		}
	}

	_constraints.initMatrices(matrixSize);

	_mesh.loadNodeProperty(_configuration.displacement       , { "X", "Y", "Z" }, { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z });
	_mesh.loadNodeProperty(_configuration.temperature        , { }              , { Property::TEMPERATURE });
	_mesh.loadNodeProperty(_configuration.obstacle           , { }              , { Property::OBSTACLE });
	_mesh.loadNodeProperty(_configuration.normal_direction   , { }              , { Property::NORMAL_DIRECTION });

	_mesh.loadProperty(_configuration.normal_presure     , { }              , { Property::PRESSURE });
	_mesh.loadProperty(_configuration.acceleration       , { "X", "Y", "Z" }, { Property::ACCELERATION_X, Property::ACCELERATION_Y, Property::ACCELERATION_Z });
	_mesh.loadProperty(_configuration.initial_temperature, { }              , { Property::INITIAL_TEMPERATURE });

	_mesh.loadMaterials(_configuration.materials, _configuration.material_set);
	_mesh.removeDuplicateRegions();
}

void Elasticity3D::saveMeshProperties(store::ResultStore &store)
{
	store.storeProperty("displacement", { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z }, store::ResultStore::ElementType::NODES);
	store.storeProperty("forces", { Property::FORCE_X, Property::FORCE_Y, Property::FORCE_Z }, store::ResultStore::ElementType::NODES);
	store.storeProperty("obstacle", { Property::OBSTACLE }, store::ResultStore::ElementType::NODES);
	store.storeProperty("normal_direction", { Property::NORMAL_DIRECTION }, store::ResultStore::ElementType::NODES);
	if (_solverConfiguration.regularization == REGULARIZATION::FIX_POINTS) {
		store::VTK::fixPoints(store.configuration(), _mesh, "fixPoints");
	}
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
		case B0_TYPE::COMBINED:
			store::VTK::mesh(store.configuration(), _mesh, "faces", store::ResultStore::ElementType::FACES);
			store::VTK::mesh(store.configuration(), _mesh, "edges", store::ResultStore::ElementType::EDGES);
			store::VTK::corners(store.configuration(), _mesh, "corners");
			break;
		case B0_TYPE::KERNELS:
			store::VTK::mesh(store.configuration(), _mesh, "faces", store::ResultStore::ElementType::FACES);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented saving properties of B0";
		}
	}
}

void Elasticity3D::saveMeshResults(store::ResultStore &store, const std::vector<std::vector<double> > &results)
{
	store.storeValues("displacement", 3, results, store::ResultStore::ElementType::NODES);
}

void Elasticity3D::assembleB1()
{
	EqualityConstraints::insertDirichletToB1(_constraints, _mesh.nodes(), pointDOFs);
	EqualityConstraints::insertElementGluingToB1(_constraints, _mesh.nodes(), pointDOFs, K);
	EqualityConstraints::insertMortarGluingToB1(_constraints, _mesh.faces(), pointDOFs);

	for (size_t i = 0; i < _mesh.evaluators().size(); i++) {
		if (_mesh.evaluators()[i]->property() == Property::OBSTACLE) {
			InequalityConstraints::insertLowerBoundToB1(_constraints, pointDOFs, { Property::OBSTACLE });
		}
	}
}

void Elasticity3D::assembleB0()
{
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		case B0_TYPE::KERNELS:
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.faces(), pointDOFs, R1);
			break;
		case B0_TYPE::COMBINED:
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.faces(), pointDOFs, R1);
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented construction of B0";
		}
	}
}

static double determinant3x3(DenseMatrix &m)
{
	const double *values = m.values();
	return fabs(
		values[0] * values[4] * values[8] +
		values[1] * values[5] * values[6] +
		values[2] * values[3] * values[7] -
		values[2] * values[4] * values[6] -
		values[1] * values[3] * values[8] -
		values[0] * values[5] * values[7]
   );
}

static void inverse(const DenseMatrix &m, DenseMatrix &inv, double det)
{
	const double *values = m.values();
	inv.resize(m.rows(), m.columns());
	double *invj = inv.values();
	double detJx = 1 / det;
	invj[0] = detJx * (values[8] * values[4] - values[7] * values[5]);
	invj[1] = detJx * (-values[8] * values[1] + values[7] * values[2]);
	invj[2] = detJx * (values[5] * values[1] - values[4] * values[2]);
	invj[3] = detJx * (-values[8] * values[3] + values[6] * values[5]);
	invj[4] = detJx * (values[8] * values[0] - values[6] * values[2]);
	invj[5] = detJx * (-values[5] * values[0] + values[3] * values[2]);
	invj[6] = detJx * (values[7] * values[3] - values[6] * values[4]);
	invj[7] = detJx * (-values[7] * values[0] + values[6] * values[1]);
	invj[8] = detJx * (values[4] * values[0] - values[3] * values[1]);
}

// B =
// dX   0   0
//  0  dY   0
//  0   0  dZ
// dY  dX   0
//  0  dZ  dY
// dZ   0  dX
static void distribute(DenseMatrix &B, DenseMatrix &dND)
{
	eslocal columns = dND.rows() * dND.columns();
	const double *dNDx = dND.values();
	const double *dNDy = dND.values() + dND.columns();
	const double *dNDz = dND.values() + 2 * dND.columns();

	double *v = B.values();

	memcpy(&v[0], dNDx,                               sizeof(double) * dND.columns());
	memcpy(&v[3 * columns + dND.columns()],     dNDx, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns + 2 * dND.columns()], dNDx, sizeof(double) * dND.columns());

	memcpy(&v[1 * columns + dND.columns()],     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[3 * columns],                     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + 2 * dND.columns()], dNDy, sizeof(double) * dND.columns());

	memcpy(&v[2 * columns + 2 * dND.columns()], dNDz, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + dND.columns()],     dNDz, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns],                     dNDz, sizeof(double) * dND.columns());
}

static void fillC(DenseMatrix &Ce, MATERIAL_MODEL model, DenseMatrix &dens, DenseMatrix &E, DenseMatrix &mi, DenseMatrix &G)
{
	switch (model) {
	case MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC: {

		double EE = E(0, 0) / ((1 + mi(0, 0)) * (1 - 2 * mi(0, 0)));

		Ce(0, 1) = Ce(0, 2) = Ce(1, 0) = Ce(1, 2) = Ce(2, 0) = Ce(2, 1) = EE * mi(0, 0);
		Ce(0, 0) = Ce(1, 1) = Ce(2, 2) = EE * (1.0 - mi(0, 0));
		Ce(3, 3) = Ce(4, 4) = Ce(5, 5) = EE * (0.5 - mi(0, 0));
		break;
	}
	case MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC:

//	       D11 = MATERIAL_Properties.D11;
//	       D12 = MATERIAL_Properties.D12;
//	       D13 = MATERIAL_Properties.D13;
//	       D14 = MATERIAL_Properties.D14;
//	       D15 = MATERIAL_Properties.D15;
//	       D16 = MATERIAL_Properties.D16;
//	       D22 = MATERIAL_Properties.D22;
//	       D23 = MATERIAL_Properties.D23;
//	       D24 = MATERIAL_Properties.D24;
//	       D25 = MATERIAL_Properties.D25;
//	       D26 = MATERIAL_Properties.D26;
//	       D33 = MATERIAL_Properties.D33;
//	       D34 = MATERIAL_Properties.D34;
//	       D35 = MATERIAL_Properties.D35;
//	       D36 = MATERIAL_Properties.D36;
//	       D44 = MATERIAL_Properties.D44;
//	       D45 = MATERIAL_Properties.D45;
//	       D46 = MATERIAL_Properties.D46;
//	       D55 = MATERIAL_Properties.D55;
//	       D56 = MATERIAL_Properties.D56;
//	       D66 = MATERIAL_Properties.D66;
//
//	       C =   [D11    D12   D13    D14   D15    D16
//	           D12    D22   D23    D24   D25    D26
//	           D13    D23   D33    D34   D35    D36
//	           D14    D24   D34    D44   D45    D46
//	           D15    D25   D35    D45   D55    D56
//	           D16    D26   D36    D46   D56    D66];

		break;

	case MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC: {

		double miXY = mi(0, 0);
		double miYZ = mi(0, 1);
		double miXZ = mi(0, 2);
		double miYX = miXY * E(0, 1) / E(0, 0);
		double miZY = miYZ * E(0, 2) / E(0, 1);
		double miZX = miXZ * E(0, 0) / E(0, 2);

		double ksi = 1 - (miXY * miYX + miYZ * miZY + miXZ * miZX) - (miXY * miYZ * miZX + miYX * miZY * miXZ);

		double dxx = E(0, 0) * (1 - miYZ * miZY) / ksi;
		double dxy = E(0, 1) * (miXY + miXZ * miZY) / ksi;
		double dxz = E(0, 2) * (miXZ + miYZ * miXY)  /ksi;
		double dyy = E(0, 1) * (1 - miXZ * miZX) / ksi;
		double dyz = E(0, 2) * (miYZ + miYX * miXZ) / ksi;
		double dzz = E(0, 2) * (1 - miYX * miXY) / ksi;


		Ce = 0;
		Ce(0, 0) = dxx; Ce(0, 1) = dxy; Ce(0, 2) = dxz;
		Ce(1, 0) = dxy; Ce(1, 1) = dyy; Ce(1, 2) = dyz;
		Ce(2, 0) = dxz; Ce(2, 1) = dyz; Ce(2, 2) = dzz;
		Ce(3, 3) = G(0, 0);
		Ce(4, 4) = G(0, 2);
		Ce(5, 5) = G(0, 1);
		break;
	}

	default:
		ESINFO(ERROR) << "Linear elasticity 3D not supports set material model";
	}
}


static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* element)
{
	DenseMatrix Ce(6, 6), XYZ(1, 3), coordinates, J, invJ, dND, B, epsilon, rhsT;
	DenseMatrix
			matDENS(element->nodes(), 1), matE(element->nodes(), 3), matMI(element->nodes(), 3), matG(element->nodes(), 3),
			matTE(element->nodes(), 3), matT(element->nodes(), 1),
			matInitT(element->nodes(), 1), inertia(element->nodes(), 3);
	DenseMatrix gpDENS(1, 1), gpE(1, 3), gpMI(1, 3), gpG(1, 3), gpTE(1, 3), gpT(1, 1), gpInitT(1, 1), gpInertia(1, 3);
	double detJ;

	const Material* material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	coordinates.resize(element->nodes(), 3);

	inertia = 0;
	for (size_t i = 0; i < element->nodes(); i++) {
		matDENS(i, 0) = material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(element->node(i));

		switch (material->getModel(PHYSICS::LINEAR_ELASTICITY_3D)) {
		case MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC:
			matE(i, 0) = matE(i, 1) = matE(i, 2) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_X)->evaluate(element->node(i));
			matMI(i, 0) = matMI(i, 1) = matMI(i, 2) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XY)->evaluate(element->node(i));
			matTE(i, 0) = matTE(i, 1) = matTE(i, 2) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_X)->evaluate(element->node(i));
			break;
		case MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC:
			ESINFO(ERROR) << "Implement ANISOTROPIC MATERIAL";
			break;
		case MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC:
			matE(i, 0) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_X)->evaluate(element->node(i));
			matE(i, 1) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_Y)->evaluate(element->node(i));
			matE(i, 2) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_Z)->evaluate(element->node(i));

			matMI(i, 0) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XY)->evaluate(element->node(i));
			matMI(i, 1) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XZ)->evaluate(element->node(i));
			matMI(i, 2) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_YZ)->evaluate(element->node(i));

			matG(i, 0) = material->get(MATERIAL_PARAMETER::SHEAR_MODULUS_XY)->evaluate(element->node(i));
			matG(i, 1) = material->get(MATERIAL_PARAMETER::SHEAR_MODULUS_XZ)->evaluate(element->node(i));
			matG(i, 2) = material->get(MATERIAL_PARAMETER::SHEAR_MODULUS_YZ)->evaluate(element->node(i));

			matTE(i, 0) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_X)->evaluate(element->node(i));
			matTE(i, 1) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_Y)->evaluate(element->node(i));
			matTE(i, 2) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_Z)->evaluate(element->node(i));
		default:
			ESINFO(ERROR) << "Linear elasticity 3D not supports set material model";
		}

		matInitT(i, 0) = element->getProperty(Property::INITIAL_TEMPERATURE, i, 0, 0);
		matT(i, 0) = element->getProperty(Property::TEMPERATURE, i, 0, matInitT(i, 0));

		inertia(i, 0) = element->sumProperty(Property::ACCELERATION_X, i, 0, 0);
		inertia(i, 1) = element->sumProperty(Property::ACCELERATION_Y, i, 0, 0);
		inertia(i, 2) = element->sumProperty(Property::ACCELERATION_Z, i, 0, 0);

		coordinates(i, 0) = mesh.coordinates()[element->node(i)].x;
		coordinates(i, 1) = mesh.coordinates()[element->node(i)].y;
		coordinates(i, 2) = mesh.coordinates()[element->node(i)].z;
	}

	eslocal Ksize = 3 * element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	std::fill(fe.begin(), fe.end(), 0);
	rhsT.resize(Ksize, 1);
	rhsT = 0;

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);
		dND.multiply(invJ, dN[gp]);
		XYZ.multiply(N[gp], coordinates);
		gpDENS.multiply(N[gp], matDENS);
		gpE.multiply(N[gp], matE);
		gpMI.multiply(N[gp], matMI);
		gpG.multiply(N[gp], matG);
		gpTE.multiply(N[gp], matTE);
		gpT.multiply(N[gp], matT);
		gpInitT.multiply(N[gp], matInitT);
		gpInertia.multiply(N[gp], inertia);

		fillC(Ce, material->getModel(PHYSICS::LINEAR_ELASTICITY_3D), gpDENS, gpE, gpMI, gpG);
		B.resize(Ce.rows(), Ksize);
		epsilon.resize(Ce.rows(), 1);

		epsilon(0, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 0);
		epsilon(1, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 1);
		epsilon(2, 0) = (gpT(0, 0) - gpInitT(0, 0)) * gpTE(0, 2);
		epsilon(3, 0) = epsilon(4, 0) = epsilon(5, 0) = 0;

		distribute(B, dND);

		Ke.multiply(B, Ce * B, detJ * weighFactor[gp], 1, true);
		rhsT.multiply(B, Ce * epsilon, detJ * weighFactor[gp], 0, true, false);
		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * N[gp](0, i % element->nodes()) * gpInertia(0, i / element->nodes());
			// fe[i] += gpDENS(0, 0) * detJ * weighFactor[gp] * N[gp](0, i % element->nodes()) * XYZ(0, i / element->nodes()) * pow(LinearElasticity2D::angularVelocity.z, 2);
			fe[i] += rhsT(i, 0);
		}
	}
}

static void processFace(std::vector<double> &fe, const espreso::Mesh &mesh, const Element* face)
{
	DenseMatrix coordinates(face->nodes(), 3), dND(1, 2), P(face->nodes(), 1), normal(3, 1), XYZ(1, 3);
	DenseMatrix gpP(1, 1), gpQ(1, 3);

	const std::vector<DenseMatrix> &dN = face->dN();
	const std::vector<DenseMatrix> &N = face->N();
	const std::vector<double> &weighFactor = face->weighFactor();

	for (size_t n = 0; n < face->nodes(); n++) {
		coordinates(n, 0) = mesh.coordinates()[face->node(n)].x;
		coordinates(n, 1) = mesh.coordinates()[face->node(n)].y;
		coordinates(n, 2) = mesh.coordinates()[face->node(n)].z;

		P(n, 0) = face->getProperty(Property::PRESSURE, n, 0, 0);
	}

	eslocal Ksize = 3 * face->nodes();
	fe.resize(Ksize);
	std::fill(fe.begin(), fe.end(), 0);

	for (size_t gp = 0; gp < face->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		face->rotateOutside(face->parentElements()[0], mesh.coordinates(), va);
		double J = va.norm();
		normal(0, 0) = va.x / va.norm();
		normal(0, 1) = va.y / va.norm();
		normal(0, 2) = va.z / va.norm();

		gpP.multiply(N[gp], P);
		gpQ.multiply(normal, gpP);

		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += J * weighFactor[gp] * N[gp](0, i % face->nodes()) * gpQ(0, i / face->nodes());
		}
	}
}

static void analyticsKernels(SparseMatrix &R1, const Coordinates &coordinates, size_t subdomain)
{
	size_t nodes = coordinates.localSize(subdomain);
	R1.rows = 3 * nodes;
	R1.cols = 6;
	R1.nnz = R1.rows * R1.cols;
	R1.type = 'G';

	R1.dense_values.reserve(R1.nnz);

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };
		kernel[c] = 1;
		for (size_t i = 0; i < nodes; i++) {
			R1.dense_values.insert(R1.dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(-p.y);
		R1.dense_values.push_back( p.x);
		R1.dense_values.push_back(   0);
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(-p.z);
		R1.dense_values.push_back(   0);
		R1.dense_values.push_back( p.x);
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(   0);
		R1.dense_values.push_back(-p.z);
		R1.dense_values.push_back( p.y);
	}
}

static void analyticsRegMat(SparseMatrix &K, SparseMatrix &RegMat, const std::vector<Element*> &fixPoints, const Coordinates &coordinates, size_t subdomain)
{
	ESTEST(MANDATORY) << "Too few FIX POINTS: " << fixPoints.size() << (fixPoints.size() > 3 ? TEST_PASSED : TEST_FAILED);

	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 6;
	Nt.cols = K.cols;
	Nt.nnz  = 9 * fixPoints.size();
	Nt.type = 'G';

	std::vector<eslocal> &ROWS = Nt.CSR_I_row_indices;
	std::vector<eslocal> &COLS = Nt.CSR_J_col_indices;
	std::vector<double>  &VALS = Nt.CSR_V_values;

	ROWS.reserve(Nt.rows + 1);
	COLS.reserve(Nt.nnz);
	VALS.reserve(Nt.nnz);

	ROWS.push_back(1);
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());

	for (size_t c = 0; c < 3; c++) {
		for (size_t i = 0; i < fixPoints.size(); i++) {
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, c) + 1);
		}
	}
	VALS.insert(VALS.end(), 3 * fixPoints.size(), 1);

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + 1);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + 1);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + 1);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 2) + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + 1);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 2) + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.y);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );
	RegMat.MatMat(Nt, 'N', N);
	RegMat.MatTranspose();
	RegMat.RemoveLower();
	RegMat.mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;

	SparseSolverCPU NtN;
	NtN.ImportMatrix(RegMat);
	RegMat.Clear();

	NtN.Factorization("Create RegMat");
	NtN.SolveMat_Sparse(Nt);
	NtN.Clear();

	RegMat.MatMat(N, 'N', Nt);
	RegMat.MatScale(K.getDiagonalMaximum());
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &RegMat, SparseMatrix &R, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);
}

void Elasticity3D::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	processElement(Ke, fe, _mesh, e);
	dofs.resize(e->nodes() * pointDOFs.size());
	for (size_t dof = 0, i = 0; dof < pointDOFs.size(); dof++) {
		for (size_t n = 0; n < e->nodes(); n++, i++) {
			dofs[i] = e->node(n) * pointDOFs.size() + dof;
		}
	}


	std::vector<Property> forces = { Property::FORCE_X, Property::FORCE_Y, Property::FORCE_Z };
	for (size_t n = 0; n < e->nodes(); n++) {
		for (size_t dof = 0; dof < pointDOFs.size(); dof++) {
			fe[n * pointDOFs.size() + dof] = e->sumProperty(forces[dof], n, 0, 0) / _mesh.nodes()[e->node(n)]->domains().size();
		}
	}
}

void Elasticity3D::makeStiffnessMatricesRegular()
{
	#pragma omp parallel for
	for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		switch (_solverConfiguration.regularization) {
		case REGULARIZATION::FIX_POINTS:
			analyticsKernels(R1[subdomain], _mesh.coordinates(), subdomain);
			analyticsRegMat(K[subdomain], RegMat[subdomain], _mesh.fixPoints(subdomain), _mesh.coordinates(), subdomain);
			K[subdomain].RemoveLower();
			RegMat[subdomain].RemoveLower();
			K[subdomain].MatAddInPlace(RegMat[subdomain], 'N', 1);
			RegMat[subdomain].ConvertToCOO(1);
			break;
		case REGULARIZATION::NULL_PIVOTS:
			K[subdomain].RemoveLower();
			algebraicKernelsAndRegularization(K[subdomain], RegMat[subdomain], R1[subdomain], subdomain);
			break;
		}
	}
}

void Elasticity3D::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(matrixSize[subdomain], matrixSize[subdomain]);
	f[subdomain].resize(matrixSize[subdomain]);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<Element*> &nodes = _mesh.nodes();

	for (eslocal e = partition[subdomain]; e < partition[subdomain + 1]; e++) {

		processElement(Ke, fe, _mesh, elements[e]);

		for (size_t nx = 0; nx < elements[e]->nodes(); nx++) {
			for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
				size_t row = nodes[elements[e]->node(nx)]->DOFIndex(subdomain, dx);
				for (size_t ny = 0; ny < elements[e]->nodes(); ny++) {
					for (size_t dy = 0; dy < pointDOFs.size(); dy++) {
						size_t column = nodes[elements[e]->node(ny)]->DOFIndex(subdomain, dy);
						_K(row, column) = Ke(dx * elements[e]->nodes() + nx, dy * elements[e]->nodes() + ny);
					}
				}
				f[subdomain][row] += fe[dx * elements[e]->nodes() + nx];
			}
		}
	}

	for (size_t i = 0; i < _mesh.faces().size(); i++) {
		if (_mesh.faces()[i]->inDomain(subdomain) && _mesh.faces()[i]->regions().size()) {
			processFace(fe, _mesh, _mesh.faces()[i]);

			for (size_t nx = 0; nx < _mesh.faces()[i]->nodes(); nx++) {
				for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
					size_t row = nodes[_mesh.faces()[i]->node(nx)]->DOFIndex(subdomain, dx);
					f[subdomain][row] += fe[dx * _mesh.faces()[i]->nodes() + nx];
				}
			}

		}
	}

	std::vector<Property> forces = { Property::FORCE_X, Property::FORCE_Y, Property::FORCE_Z };
	for (size_t n = 0; n < _mesh.coordinates().localSize(subdomain); n++) {
		Element *node = _mesh.nodes()[_mesh.coordinates().clusterIndex(n, subdomain)];
		for (size_t dof = 0; dof < pointDOFs.size(); dof++) {
			f[subdomain][node->DOFIndex(subdomain, dof)] += node->sumProperty(forces[dof], 0, 0, 0) / node->numberOfGlobalDomainsWithDOF(dof);
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
}

static void postProcessElement(std::vector<double> &stress, std::vector<double> &principleStress, std::vector<double> &HMH, DenseMatrix &solution, const Element* element, const Mesh &mesh, const LinearElasticity3DConfiguration &configuration)
{
	DenseMatrix Ce(6, 6), XYZ(1, 3), coordinates, J, invJ, dND, B, epsilon, rhsT;
	DenseMatrix
			matDENS(element->nodes(), 1), matE(element->nodes(), 3), matMI(element->nodes(), 3), matG(element->nodes(), 3),
			matTE(element->nodes(), 3), matT(element->nodes(), 1),
			matInitT(element->nodes(), 1), inertia(element->nodes(), 3);
	DenseMatrix gpDENS(1, 1), gpE(1, 3), gpMI(1, 3), gpG(1, 3), gpTE(1, 3), gpT(1, 1), gpInitT(1, 1), gpInertia(1, 3);
	DenseMatrix matStress(6, 1);
	double detJ;

	const Material* material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();

	coordinates.resize(element->nodes(), 3);

	inertia = 0;
	for (size_t i = 0; i < element->nodes(); i++) {
		matDENS(i, 0) = material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(element->node(i));

		switch (material->getModel(PHYSICS::LINEAR_ELASTICITY_3D)) {
		case MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC:
			matE(i, 0) = matE(i, 1) = matE(i, 2) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_X)->evaluate(element->node(i));
			matMI(i, 0) = matMI(i, 1) = matMI(i, 2) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XY)->evaluate(element->node(i));
			matTE(i, 0) = matTE(i, 1) = matTE(i, 2) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_X)->evaluate(element->node(i));
			break;
		case MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC:
			ESINFO(ERROR) << "Implement ANISOTROPIC MATERIAL";
			break;
		case MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC:
			matE(i, 0) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_X)->evaluate(element->node(i));
			matE(i, 1) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_Y)->evaluate(element->node(i));
			matE(i, 2) = material->get(MATERIAL_PARAMETER::YOUNG_MODULUS_Z)->evaluate(element->node(i));

			matMI(i, 0) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XY)->evaluate(element->node(i));
			matMI(i, 1) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_XZ)->evaluate(element->node(i));
			matMI(i, 2) = material->get(MATERIAL_PARAMETER::POISSON_RATIO_YZ)->evaluate(element->node(i));

			matG(i, 0) = material->get(MATERIAL_PARAMETER::SHEAR_MODULUS_XY)->evaluate(element->node(i));
			matG(i, 1) = material->get(MATERIAL_PARAMETER::SHEAR_MODULUS_XZ)->evaluate(element->node(i));
			matG(i, 2) = material->get(MATERIAL_PARAMETER::SHEAR_MODULUS_YZ)->evaluate(element->node(i));

			matTE(i, 0) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_X)->evaluate(element->node(i));
			matTE(i, 1) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_Y)->evaluate(element->node(i));
			matTE(i, 2) = material->get(MATERIAL_PARAMETER::THERMAL_EXPANSION_Z)->evaluate(element->node(i));
		default:
			ESINFO(ERROR) << "Linear elasticity 3D not supports set material model";
		}

		matInitT(i, 0) = element->getProperty(Property::INITIAL_TEMPERATURE, i, 0, 0);
		matT(i, 0) = element->getProperty(Property::TEMPERATURE, i, 0, matInitT(i, 0));

		coordinates(i, 0) = mesh.coordinates()[element->node(i)].x;
		coordinates(i, 1) = mesh.coordinates()[element->node(i)].y;
		coordinates(i, 2) = mesh.coordinates()[element->node(i)].z;
	}

	eslocal Ksize = 3 * element->nodes();

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);
		dND.multiply(invJ, dN[gp]);
		gpDENS.multiply(N[gp], matDENS);
		gpE.multiply(N[gp], matE);
		gpMI.multiply(N[gp], matMI);
		gpG.multiply(N[gp], matG);

		fillC(Ce, material->getModel(PHYSICS::LINEAR_ELASTICITY_3D), gpDENS, gpE, gpMI, gpG);
		B.resize(Ce.rows(), Ksize);
		distribute(B, dND);

		matStress.multiply(Ce, B * solution, 1, 1);
	}
	for (size_t i = 0; i < 6; i++) {
		stress.push_back(matStress(i, 0) / element->gaussePoints());
	}

	std::vector<double> tensor = {
			matStress(0, 0) / element->gaussePoints(),
			matStress(3, 0) / element->gaussePoints(),
			matStress(5, 0) / element->gaussePoints(),
			matStress(1, 0) / element->gaussePoints(),
			matStress(4, 0) / element->gaussePoints(),
			matStress(2, 0) / element->gaussePoints()
	};
	principleStress.insert(principleStress.end(), 3, 0);
	LAPACKE_dspev(LAPACK_ROW_MAJOR, 'N', 'U', 3, tensor.data(), principleStress.data() + principleStress.size() - 3, NULL, 3);

	auto hmh = [] (double *data) {
		return sqrt((pow(data[0] - data[1], 2) + pow(data[0] - data[2], 2) + pow(data[1] - data[2], 2)) / 2);
	};

	HMH.push_back(hmh(principleStress.data() + principleStress.size() - 3));
}

void Elasticity3D::postProcess(store::ResultStore &store, const std::vector<std::vector<double> > &solution)
{
	if (!_configuration.post_process) {
		return;
	}
	std::vector<std::vector<double> > stress(_mesh.parts());
	std::vector<std::vector<double> > principleStress(_mesh.parts());
	std::vector<std::vector<double> > HMH(_mesh.parts());
	DenseMatrix eSolution;

	for (size_t p = 0; p < _mesh.parts(); p++) {
		stress[p].reserve(6 * matrixSize[p]);
		principleStress[p].reserve(3 * matrixSize[p]);
		HMH[p].reserve(1 * matrixSize[p]);
		for (eslocal e = _mesh.getPartition()[p]; e < _mesh.getPartition()[p + 1]; e++) {
			eSolution.resize(pointDOFs.size() * _mesh.elements()[e]->nodes(), 1);
			for (size_t dof = 0, i = 0; dof < pointDOFs.size(); dof++) {
				for (size_t n = 0; n < _mesh.elements()[e]->nodes(); n++, i++) {
					eSolution(i, 0) = solution[p][_mesh.nodes()[_mesh.elements()[e]->node(n)]->DOFIndex(p, dof)];
				}
			}
			postProcessElement(stress[p], principleStress[p], HMH[p], eSolution, _mesh.elements()[e], _mesh, _configuration);
		}
	}

	store.storeValues("stress", 6, stress, store::ResultStore::ElementType::ELEMENTS);
	store.storeValues("principle_stress", 3, principleStress, store::ResultStore::ElementType::ELEMENTS);
	store.storeValues("HMH", 1, HMH, store::ResultStore::ElementType::ELEMENTS);
}

}


