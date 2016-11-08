
#include "assembler.h"

using namespace espreso;

std::vector<Property> LinearElasticity3D::elementDOFs;
std::vector<Property> LinearElasticity3D::faceDOFs;
std::vector<Property> LinearElasticity3D::edgeDOFs;
std::vector<Property> LinearElasticity3D::pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };
std::vector<Property> LinearElasticity3D::midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };

void LinearElasticity3D::prepareMeshStructures()
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

	if (config::solver::REGULARIZATION == config::solver::REGULARIZATIONalternative::FIX_POINTS) {
		_mesh.computeFixPoints(config::mesh::FIX_POINTS);
	}

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			_mesh.computeVolumeCorners(config::mesh::CORNERS, config::mesh::VERTEX_CORNERS, config::mesh::EDGE_CORNERS, config::mesh::FACE_CORNERS);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			_mesh.computeFacesSharedByDomains();
			break;
		case config::solver::B0_TYPEalternative::COMBINED:
			_mesh.computeFacesSharedByDomains();
			if (!_mesh.corners().size()) {
				_mesh.computeEdgesOnBordersOfFacesSharedByDomains();
				_mesh.computeCornersOnEdges(config::mesh::VERTEX_CORNERS ? config::mesh::CORNERS : 0);
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented B0";
		}
	}
}

void LinearElasticity3D::saveMeshProperties(output::Store &store)
{
	store.storeProperty("displacement", { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z }, output::Store::ElementType::NODES);
	store.storeProperty("forces", { Property::FORCE_X, Property::FORCE_Y, Property::FORCE_Z }, output::Store::ElementType::NODES);
	store.storeProperty("obstacle", { Property::OBSTACLE }, output::Store::ElementType::NODES);
	store.storeProperty("normal_direction", { Property::NORMAL_DIRECTION }, output::Store::ElementType::NODES);
	if (config::solver::REGULARIZATION == config::solver::REGULARIZATIONalternative::FIX_POINTS) {
		output::VTK::fixPoints(_mesh, "fixPoints", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
		case config::solver::B0_TYPEalternative::COMBINED:
			output::VTK::mesh(_mesh, "faces", output::Store::ElementType::FACES, config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
			output::VTK::mesh(_mesh, "edges", output::Store::ElementType::EDGES, config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
			output::VTK::corners(_mesh, "corners", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			output::VTK::mesh(_mesh, "faces", output::Store::ElementType::FACES, config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented saving properties of B0";
		}
	}
}

void LinearElasticity3D::saveMeshResults(output::Store &store, const std::vector<std::vector<double> > &results)
{
	store.storeValues("displacement", 3, results, output::Store::ElementType::NODES);
	store.finalize();
}

void LinearElasticity3D::assembleGluingMatrices()
{
	_constraints.initMatrices(matrixSize);

	EqualityConstraints::insertDirichletToB1(_constraints, _mesh.nodes(), pointDOFs);
	EqualityConstraints::insertElementGluingToB1(_constraints, _mesh.nodes(), pointDOFs, K);
	EqualityConstraints::insertMortarGluingToB1(_constraints, _mesh.faces(), pointDOFs);

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.faces(), pointDOFs, R1);
			break;
		case config::solver::B0_TYPEalternative::COMBINED:
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.faces(), pointDOFs, R1);
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented construction of B0";
		}
	}

	for (size_t i = 0; i < _mesh.evaluators().size(); i++) {
		if (_mesh.evaluators()[i]->property() == Property::OBSTACLE) {
			InequalityConstraints::insertLowerBoundToB1(_constraints, _mesh.nodes(), pointDOFs, { Property::OBSTACLE });
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

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* element)
{
	DenseMatrix Ce(6, 6), coordinates, J, invJ, dND, B;
	std::vector<double> inertia(3, 0);
	double detJ;

	const Material &material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	// TODO: set the omega from example
	Point omega(50, 50, 0);

	double ex = material.youngModulusX(0);
	double mi = material.poissonRatioXY(0);
	double E = ex / ((1 + mi) * (1 - 2 * mi));
	Ce(0, 1) = Ce(0, 2) = Ce(1, 0) = Ce(1, 2) = Ce(2, 0) = Ce(2, 1) = E * mi;
	Ce(0, 0) = Ce(1, 1) = Ce(2, 2) = E * (1.0 - mi);
	Ce(3, 3) = Ce(4, 4) = Ce(5, 5) = E * (0.5 - mi);

	inertia[0] = inertia[1] = 0; // inertia[2] = 0;
	inertia[2] = 9.8066 * material.density(0);

	coordinates.resize(element->nodes(), 3);

	Point mid;
	for (size_t i = 0; i < element->nodes(); i++) {
		coordinates(i, 0) = mesh.coordinates()[element->node(i)].x;
		coordinates(i, 1) = mesh.coordinates()[element->node(i)].y;
		coordinates(i, 2) = mesh.coordinates()[element->node(i)].z;
		mid += mesh.coordinates()[element->node(i)];
	}
	mid /= element->nodes();

	eslocal Ksize = 3 * element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	double rotation[3] = { mid.x * omega.x * omega.x, mid.y * omega.y * omega.y, mid.z * omega.z * omega.z };

	for (eslocal gp = 0; gp < element->gaussePoints(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);
		B.resize(Ce.rows(), Ksize);
		distribute(B, dND);
		Ke.multiply(B, Ce * B, detJ * weighFactor[gp], 1, true);

		for (eslocal i = 0; i < Ksize; i++) {
			// TODO: set rotation from example
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i % element->nodes()) * inertia[i / element->nodes()];
			//fe[i] += detJ * weighFactor[gp] * N[gp](0, i % e->size()) * 7850 * rotation[i / e->size()];
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
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, c) + IJVMatrixIndexing);
		}
	}
	VALS.insert(VALS.end(), 3 * fixPoints.size(), 1);

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + IJVMatrixIndexing);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + IJVMatrixIndexing);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + IJVMatrixIndexing);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 2) + IJVMatrixIndexing);
		VALS.push_back(-p.z);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + IJVMatrixIndexing);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 2) + IJVMatrixIndexing);
		VALS.push_back(-p.z);
		VALS.push_back( p.y);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );
	RegMat.MatMat(Nt, 'N', N);
	RegMat.MatTranspose();
	RegMat.RemoveLower();
	RegMat.mtype = SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;

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

void LinearElasticity3D::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
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
			if (_mesh.nodes()[e->node(n)]->settings().isSet(forces[dof])) {
				fe[n * pointDOFs.size() + dof] = _mesh.nodes()[e->node(n)]->settings(forces[dof]).back()->evaluate(e->node(n)) / _mesh.nodes()[e->node(n)]->domains().size();
			}
		}
	}
}

void LinearElasticity3D::makeStiffnessMatricesRegular()
{
	cilk_for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		switch (config::solver::REGULARIZATION) {
		case config::solver::REGULARIZATIONalternative::FIX_POINTS:
			if (config::assembler::DOFS_ORDER == config::assembler::DOFS_ORDERalternative::GROUP_DOFS) {
				ESINFO(GLOBAL_ERROR) << "Implement regularization for GROUP_DOFS alternative";
			}

			analyticsKernels(R1[subdomain], _mesh.coordinates(), subdomain);
			analyticsRegMat(K[subdomain], RegMat[subdomain], _mesh.fixPoints(subdomain), _mesh.coordinates(), subdomain);
			K[subdomain].RemoveLower();
			RegMat[subdomain].RemoveLower();
			K[subdomain].MatAddInPlace(RegMat[subdomain], 'N', 1);
			RegMat[subdomain].ConvertToCOO(1);
			break;
		case config::solver::REGULARIZATIONalternative::NULL_PIVOTS:
			K[subdomain].RemoveLower();
			algebraicKernelsAndRegularization(K[subdomain], RegMat[subdomain], R1[subdomain], subdomain);
			break;
		}
	}
}

void LinearElasticity3D::composeSubdomain(size_t subdomain)
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

	std::vector<Property> forces = { Property::FORCE_X, Property::FORCE_Y, Property::FORCE_Z };
	for (size_t n = 0; n < _mesh.coordinates().localSize(subdomain); n++) {
		Element *node = _mesh.nodes()[_mesh.coordinates().clusterIndex(n, subdomain)];
		for (size_t dof = 0; dof < pointDOFs.size(); dof++) {
			if (node->settings().isSet(forces[dof])) {
				f[subdomain][node->DOFIndex(subdomain, dof)] = node->settings(forces[dof]).back()->evaluate(node->node(0)) / node->domains().size();
			}
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
}




