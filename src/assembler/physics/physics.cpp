
#include "physics.h"

#include "../../basis/containers/serializededata.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/parser.h"
#include "../../basis/utilities/communication.h"
#include "../instance.h"
#include "../step.h"

#include "../../mesh/mesh.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/store/nodestore.h"
#include "../../mesh/store/elementstore.h"
#include "../../mesh/store/boundaryregionstore.h"
#include "../../mesh/store/elementsregionstore.h"

#include "../constraints/equalityconstraints.h"

#include "../../solver/generic/SparseMatrix.h"
#include "../../basis/matrices/sparseVVPMatrix.h"
#include "../../basis/matrices/denseMatrix.h"
#include "../../basis/matrices/sparseCSRMatrix.h"
#include "../../config/ecf/solver/feti.h"
#include "../../config/ecf/physics/physics.h"

#ifdef BEM4I
#include "esbem.h"
#endif

using namespace espreso;

Physics::Physics()
: _name(""), _mesh(NULL), _instance(NULL), _step(NULL), _equalityConstraints(NULL), _configuration(NULL), _DOFs(0)
{

}

Physics::Physics(const std::string &name, Mesh *mesh, Instance *instance, Step *step, const PhysicsConfiguration *configuration, int DOFs)
: _name(name), _mesh(mesh), _instance(instance), _step(step), _equalityConstraints(NULL), _configuration(configuration), _DOFs(DOFs) // initialized in a particular physics
{
	std::vector<int> BEMRegions(_mesh->elements->regionMaskSize);
	for (auto it = configuration->discretization.begin(); it != configuration->discretization.end(); ++it) {
		if (it->second == DISCRETIZATION::BEM) {
			for (size_t r = 0; r < _mesh->elementsRegions.size(); r++) {
				if (StringCompare::caseInsensitiveEq(it->first, _mesh->elementsRegions[r]->name)) {
					BEMRegions[r / (8 * sizeof(int))] |= 1 << (r % (8 * sizeof(int)));
				}
			}
		}
	}

	_hasBEM = false;
	_BEMDomain.resize(_mesh->elements->ndomains);
	for (eslocal d = 0; d < _mesh->elements->ndomains; d++) {
		for (int i = 0; i < _mesh->elements->regionMaskSize; i++) {
			if (_mesh->elements->regions->datatarray()[_mesh->elements->elementsDistribution[d] * _mesh->elements->regionMaskSize + i] & BEMRegions[i]) {
				_BEMDomain[d] = 1;
				_hasBEM = true;
			}
		}
	}

	_BEMData.resize(mesh->elements->ndomains, NULL);
}

Physics::~Physics()
{
	if (_equalityConstraints != NULL) {
		delete _equalityConstraints;
	}
#ifdef BEM4I
	for (size_t i = 0; i < _BEMData.size(); i++) {
		if (_BEMData[i] != NULL) {
			bem4i::deleteBem4iData(_BEMData[i]);
		}
	}
#endif
}

void Physics::printInvalidElement(eslocal eindex) const
{
	auto nodes = _mesh->elements->nodes->begin() + eindex;

	std::ofstream os("invalidElement.vtk");
	os << "# vtk DataFile Version 2.0\n";
	os << "INVALID ELEMENT\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n";
	os << "\n";
	os << "POINTS " << nodes->size() << " float\n";
	for (auto n = nodes->begin(); n != nodes->end(); ++n) {
		os << _mesh->nodes->coordinates->datatarray()[*n].x << " " << _mesh->nodes->coordinates->datatarray()[*n].y << " " << _mesh->nodes->coordinates->datatarray()[*n].z << "\n";
	}
	os << "\n";
	os << "CELLS 1 " << nodes->size() + 1 << "\n";
	os << nodes->size();
	for (auto n = nodes->begin(); n != nodes->end(); ++n) {
		os << " " << n - nodes->begin();
	}
	os << "\n";
	os << "CeLL_TYPES 1\n";
	switch (_mesh->elements->epointers->datatarray()[eindex]->code) {
	case Element::CODE::SQUARE4:
		os << "9\n";
		break;
	case Element::CODE::SQUARE8:
		os << "23\n";
		break;
	case Element::CODE::TRIANGLE3:
		os << "5\n";
		break;
	case Element::CODE::TRIANGLE6:
		os << "22\n";
		break;
	case Element::CODE::TETRA4:
		os << "10\n";
		break;
	case Element::CODE::TETRA10:
		os << "24\n";
		break;
	case Element::CODE::PYRAMID5:
		os << "14\n";
		break;
	case Element::CODE::PYRAMID13:
		os << "27\n";
		break;
	case Element::CODE::PRISMA6:
		os << "13\n";
		break;
	case Element::CODE::PRISMA15:
		os << "26\n";
		break;
	case Element::CODE::HEXA8:
		os << "12\n";
		break;
	case Element::CODE::HEXA20:
		os << "25\n";
		break;
	default:
		break;
	}
}

void Physics::updateMatrix(Matrices matrix)
{
	#pragma omp parallel for
	for  (size_t d = 0; d < _instance->domains; d++) {

		updateMatrix(matrix, d);

		switch (_instance->K[d].mtype) {
		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
			_instance->K[d].RemoveLower();
			_instance->M[d].RemoveLower();
			break;
		case MatrixType::REAL_UNSYMMETRIC:
			break;
		}

		ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	ESINFO(PROGRESS3);
}

void Physics::updateMatrix(Matrices matrices, size_t domain)
{
	SparseVVPMatrix<eslocal> _K, _M;

	if (matrices & Matrices::f) {
		_instance->f[domain].clear();
		_instance->f[domain].resize(_instance->domainDOFCount[domain]);
	}

	if (_BEMDomain[domain]) {
		if (matrices & Matrices::M) {
			ESINFO(ERROR) << "BEM not support computation of matrix M.";
		}
		if (matrices & Matrices::R) {
			ESINFO(ERROR) << "BEM not support computation of matrix R.";
		}
		processBEM(domain, matrices);
		_instance->K[domain].ConvertDenseToCSR(1);
		assembleBoundaryConditions(_K, _M, matrices & Matrices::f, domain);
	} else {
		if (matrices & Matrices::K) {
			_K.resize(_instance->domainDOFCount[domain], _instance->domainDOFCount[domain]);
		}
		if (matrices & Matrices::M) {
			_M.resize(_instance->domainDOFCount[domain], _instance->domainDOFCount[domain]);
		}
		if (matrices & Matrices::R) {
			_instance->R[domain].clear();
			_instance->R[domain].resize(_instance->domainDOFCount[domain]);
		}

		std::vector<eslocal> DOFs;
		DenseMatrix Ke, Me, Re, fe;

		auto nodes = _mesh->elements->nodes->cbegin() + _mesh->elements->elementsDistribution[domain];
		for (eslocal e = _mesh->elements->elementsDistribution[domain]; e < (eslocal)_mesh->elements->elementsDistribution[domain + 1]; ++e, ++nodes) {
			processElement(domain, matrices, e, Ke, Me, Re, fe);
			fillDOFsIndices(*nodes, domain, DOFs);
			insertElementToDomain(_K, _M, DOFs, Ke, Me, Re, fe, domain, false);
		}
		assembleBoundaryConditions(_K, _M, matrices, domain);
		// TODO: make it direct
		if (matrices & Matrices::K) {
			SparseCSRMatrix<eslocal> csrK = _K;
			_instance->K[domain] = csrK;
			_instance->K[domain].mtype = getMatrixType(domain);
		}
		if (matrices & Matrices::M) {
			SparseCSRMatrix<eslocal> csrM = _M;
			_instance->M[domain] = csrM;
			_instance->M[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
		}
	}
}

void Physics::assembleBoundaryConditions(SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, Matrices matrices, size_t domain)
{
	DenseMatrix Ke, fe, Me(0, 0), Re(0, 0);
	std::vector<eslocal> DOFs;

	for (size_t r = 0; r < _mesh->boundaryRegions.size(); r++) {
		if (_mesh->boundaryRegions[r]->dimension == 2) {
			if (_mesh->boundaryRegions[r]->eintervalsDistribution[domain] < _mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
				eslocal begin = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
				eslocal end = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
				auto nodes = _mesh->boundaryRegions[r]->elements->cbegin() + begin;
				for (eslocal i = begin; i < end; ++i, ++nodes) {
					processFace(domain, _mesh->boundaryRegions[r], matrices, i, Ke, Me, Re, fe);
					fillDOFsIndices(*nodes, domain, DOFs);
					insertElementToDomain(K, M, DOFs, Ke, Me, Re, fe, domain, true);
				}
			}
		}
		if (_mesh->boundaryRegions[r]->dimension == 1) {
			if (_mesh->boundaryRegions[r]->eintervalsDistribution[domain] < _mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
				eslocal begin = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
				eslocal end = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
				auto nodes = _mesh->boundaryRegions[r]->elements->cbegin() + begin;
				for (eslocal i = begin; i < end; ++i, ++nodes) {
					processEdge(domain, _mesh->boundaryRegions[r], matrices, i, Ke, Me, Re, fe);
					fillDOFsIndices(*nodes, domain, DOFs);
					insertElementToDomain(K, M, DOFs, Ke, Me, Re, fe, domain, true);
				}
			}
		}
		// TODO: process NODE
//		if (_mesh->boundaryRegions[r]->dimension == 0) {
//			if (_mesh->boundaryRegions[r]->eintervalsDistribution[domain] < _mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
//				eslocal begin = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
//				eslocal end = _mesh->boundaryRegions[r]->eintervals[_mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
//				auto nodes = _mesh->boundaryRegions[r]->elements->cbegin() + begin;
//				for (eslocal i = begin; i < end; ++i, ++nodes) {
//					processEdge(domain, _mesh->boundaryRegions[r], matrices, i, Ke, Me, Re, fe);
//					fillDOFsIndices(*nodes, domain, DOFs);
//					insertElementToDomain(K, M, DOFs, Ke, Me, Re, fe, domain, true);
//				}
//			}
//		}
	}
}

/**
 *
 * The method assumed that element matrix is composed in the following order:
 * x1, x2, x3, ..., y1, y2, y3, ..., z1, z2, z3,...
 *
 */
void Physics::fillDOFsIndices(edata<const eslocal> &nodes, eslocal domain, std::vector<eslocal> &DOFs) const
{
	// TODO: improve performance
	const std::vector<DomainInterval> &intervals = _mesh->nodes->dintervals[domain];
	DOFs.resize(_DOFs * nodes.size());
	size_t i = 0;
	for (size_t dof = 0; dof < _DOFs; dof++) {
		for (auto n = nodes.begin(); n != nodes.end(); n++) {
			auto it = std::lower_bound(intervals.begin(), intervals.end(), *n, [] (const DomainInterval &interval, eslocal node) { return interval.end <= node; });
			DOFs[i++] = _DOFs * (it->DOFOffset + *n - it->begin) + dof;
		}
	}
}

void Physics::insertElementToDomain(
		SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M,
		const std::vector<eslocal> &DOFs,
		const DenseMatrix &Ke, const DenseMatrix &Me, const DenseMatrix &Re, const DenseMatrix &fe,
		size_t domain, bool isBOundaryCondition)
{
	double RHSreduction = _step->internalForceReduction;
	double Kreduction = isBOundaryCondition ? RHSreduction : 1;

	if (Ke.rows() == DOFs.size() && Ke.columns() == DOFs.size()) {
		for (size_t r = 0; r < DOFs.size(); r++) {
			for (size_t c = 0; c < DOFs.size(); c++) {
				K(DOFs[r], DOFs[c]) = Kreduction * Ke(r, c);
			}
		}
	} else {
		if (Ke.rows() != 0 || Ke.columns() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while assembling stiffness matrix K(" << Ke.rows() << "," << Ke.columns() << ").";
		}
	}

	if (Me.rows() && Me.columns() && DOFs.size() % Me.rows() == 0 && DOFs.size() % Me.columns() == 0) {
		size_t multiplicity = DOFs.size() / Me.rows();
		for (size_t m = 0; m < multiplicity; m++) {
			for (size_t r = 0; r < Me.rows(); r++) {
				for (size_t c = 0; c < Me.columns(); c++) {
					M(DOFs[r * multiplicity + m], DOFs[c * multiplicity + m]) = Me(r, c);
				}
			}
		}
	} else {
		if (Me.rows() != 0 || Me.columns() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while assembling mass matrix M(" << Me.rows() << "," << Me.columns() << ").";
		}
	}


	if (Re.rows() == DOFs.size() && Re.columns() == 1) {
		for (size_t r = 0; r < DOFs.size(); r++) {
			_instance->R[domain][DOFs[r]] += Re(r, 0);
		}
	} else {
		if (Re.rows() != 0 || Re.columns() != 0) {
			std::cout << DOFs;
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while assembling matrix R(" << Re.rows() << "," << Re.columns() << ") with residual forces.";
		}
	}

	if (fe.rows() == DOFs.size() && fe.columns() == 1) {
		for (size_t r = 0; r < DOFs.size(); r++) {
			_instance->f[domain][DOFs[r]] += RHSreduction * fe(r, 0);
		}
	} else {
		if (fe.rows() != 0 || fe.columns() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: something wrong happens while assembling right-hand side vector f(" << fe.rows() << "," << fe.columns() << ").";
		}
	}
}

void Physics::makeStiffnessMatricesRegular(FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster)
{
	#pragma omp parallel for
	for (size_t d = 0; d < _instance->domains; d++) {
		makeStiffnessMatrixRegular(regularization, scSize, d, ortogonalCluster);
		ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	ESINFO(PROGRESS3);
}

void Physics::makeStiffnessMatrixRegular(FETI_REGULARIZATION regularization, size_t scSize, size_t domain, bool ortogonalCluster)
{
	switch (regularization) {

	case FETI_REGULARIZATION::ANALYTIC:
		analyticRegularization(domain, ortogonalCluster);
		_instance->RegMat[domain].RemoveLower();
		_instance->K[domain].MatAddInPlace(_instance->RegMat[domain], 'N', 1);
		_instance->RegMat[domain].ConvertToCOO(1);
		break;

	case FETI_REGULARIZATION::ALGEBRAIC:
		switch (_instance->K[domain].mtype) {
			double norm;
			eslocal defect;

		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			_instance->K[domain].get_kernel_from_K(_instance->K[domain], _instance->RegMat[domain], _instance->N1[domain], norm, defect, domain, scSize);
			break;

		case MatrixType::REAL_UNSYMMETRIC:
			_instance->K[domain].get_kernels_from_nonsym_K(_instance->K[domain], _instance->RegMat[domain], _instance->N1[domain], _instance->N2[domain], norm, defect, domain, scSize);
			break;

		default:
			ESINFO(ERROR) << "Unknown matrix type for regularization.";
		}
		break;
	}
}

double Physics::sumSquares(const std::vector<std::vector<double> > &data, SumRestriction restriction) const
{
	ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement sumSquares.";
	return 0;
}

void Physics::assembleB1(bool withRedundantMultipliers, bool withGluing, bool withScaling)
{
	_equalityConstraints->B1DirichletInsert(*_step);
	if (withGluing) {
		_equalityConstraints->B1GlueElements();
	}
}

void Physics::updateDirichletInB1(bool withRedundantMultipliers)
{
	_equalityConstraints->B1DirichletUpdate(*_step);
}

void Physics::assembleB0FromCorners()
{
	_equalityConstraints->B0Corners();
}

void Physics::assembleB0FromKernels(const std::vector<SparseMatrix> &kernels)
{
	_equalityConstraints->B0Kernels(kernels);
}

void Physics::smoothstep(double &smoothStep, double &derivation, double edge0, double edge1, double value, size_t order)
{
	value = std::max(0.0, std::min((value - edge0) / (edge1 - edge0), 1.0));

	switch (order) {
	case 0:
		smoothStep = value;
		if (value == 0 || value == 1) {
			derivation = 0;
		} else {
			derivation = 1;
		}
		break;
	case 1:
		smoothStep = -2 * pow(value, 3) + 3 * pow(value, 2);
		derivation = -6 * pow(value, 2) + 6 * value;
		break;
	case 2:
		smoothStep = 6 * pow(value, 5) - 15 * pow(value, 4) + 10 * pow(value, 3);
		derivation = 30 * pow(value, 4) - 60*  pow(value, 3) + 30 * pow(value, 2);
		break;
	case 3:
		smoothStep = -20 * pow(value, 7) + 70 * pow(value, 6) - 84 * pow(value, 5) + 35 * pow(value, 4);
		derivation = -140 * pow(value, 6) + 420 * pow(value, 5) - 420 * pow(value, 4) + 140 * pow(value, 3);
		break;
	case 4:
		smoothStep = 70 * pow(value, 9) - 315 * pow(value, 8) + 540 * pow(value, 7) - 420 * pow(value, 6) + 126 * pow(value, 5);
		derivation = 630 * pow(value, 8) - 2520 * pow(value, 7) + 3780 * pow(value, 6) - 2520 * pow(value, 5) + 630 * pow(value, 4);
		break;
	case 5:
		smoothStep = -252 * pow(value, 11) + 1386 * pow(value, 10) - 3080 * pow(value, 9) + 3465 * pow(value, 8) - 1980 * pow(value, 7) + 462 * pow(value, 6);
		derivation = -2772 * pow(value, 10) + 13860 * pow(value, 9) - 27720 * pow(value, 8) + 27720 * pow(value, 7) - 13860 * pow(value, 6) + 2772 * pow(value, 5);
		break;
	case 6:
		smoothStep = 924 * pow(value, 13) - 6006 * pow(value, 12) + 16380 * pow(value, 11) - 24024 * pow(value, 10) + 20020 * pow(value, 9) - 9009 * pow(value, 8) + 1716 * pow(value, 7);
		derivation = 12012 * pow(value, 12) - 72072 * pow(value, 11) + 180180 * pow(value, 10) - 240240 * pow(value, 9) + 180180 * pow(value, 8) - 72072 * pow(value, 7) + 12012 * pow(value, 6);
		break;
	}

	derivation /= edge1 - edge0;
}


