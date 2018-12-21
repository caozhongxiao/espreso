
#include "structuralmechanics3d.fetiprovider.h"

#include "../../../dataholder.h"

#include "../../../../globals/run.h"

#include "../../../../basis/containers/serializededata.h"
#include "../../../../basis/matrices/matrixtype.h"
#include "../../../../basis/logging/logging.h"
#include "../../../../config/ecf/physics/structuralmechanics.h"

#include "../../../../mesh/mesh.h"
#include "../../../../mesh/store/elementstore.h"
#include "../../../../mesh/store/nodestore.h"
#include "../../../../mesh/store/fetidatastore.h"

#include "../../../../solver/generic/SparseMatrix.h"
#include "../../../../solver/specific/sparsesolvers.h"

using namespace espreso;

StructuralMechanics3DFETIProvider::StructuralMechanics3DFETIProvider(StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanicsFETIProvider(configuration)
{

}

void StructuralMechanics3DFETIProvider::analyticRegularization(esint domain, bool ortogonalCluster)
{
	if (run::data->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC";
	}

	Point center = _dCenter[domain], norm = _dNorm[domain];
	double r44 = _dr44[domain], r45 = _dr45[domain], r46 = _dr46[domain], r55 = _dr55[domain], r56 = _dr56[domain];
	size_t np = _dNp[domain];

	if (ortogonalCluster) {
		size_t cluster = run::mesh->elements->clusters[domain];
		center = _cCenter[cluster], norm = _cNorm[cluster];
		r44 = _cr44[cluster], r45 = _cr45[cluster], r46 = _cr46[cluster], r55 = _cr55[cluster], r56 = _cr56[cluster];
		np = _cNp[cluster];
	} else {
		center = _dCenter[domain], norm = _dNorm[domain];
		r44 = _dr44[domain], r45 = _dr45[domain], r46 = _dr46[domain], r55 = _dr55[domain], r56 = _dr56[domain];
		np = _dNp[domain];
	}

	run::data->N1[domain].rows = run::data->K[domain].rows;
	run::data->N1[domain].cols = 6;
	run::data->N1[domain].nnz = run::data->N1[domain].rows * run::data->N1[domain].cols;
	run::data->N1[domain].type = 'G';

	run::data->N1[domain].dense_values.reserve(run::data->N1[domain].nnz);

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };
		kernel[c] = 1 / std::sqrt(np);
		for (size_t i = 0; i < run::data->K[domain].rows / 3; i++) {
			run::data->N1[domain].dense_values.insert(run::data->N1[domain].dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < run::mesh->nodes->dintervals[domain].size(); i++) {
		for (esint n = run::mesh->nodes->dintervals[domain][i].begin; n < run::mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = run::mesh->nodes->coordinates->datatarray()[n] - center;
			run::data->N1[domain].dense_values.push_back(-p.y / norm.x);
			run::data->N1[domain].dense_values.push_back( p.x / norm.x);
			run::data->N1[domain].dense_values.push_back(             0);
		}
	}

	for (size_t i = 0; i < run::mesh->nodes->dintervals[domain].size(); i++) {
		for (esint n = run::mesh->nodes->dintervals[domain][i].begin; n < run::mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = run::mesh->nodes->coordinates->datatarray()[n] - center;
			run::data->N1[domain].dense_values.push_back((-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y);
			run::data->N1[domain].dense_values.push_back((   0 - r45 / r44 * ( p.x / norm.x)) / norm.y);
			run::data->N1[domain].dense_values.push_back(( p.x - r45 / r44 * (   0 / norm.x)) / norm.y);
		}
	}

	for (size_t i = 0; i < run::mesh->nodes->dintervals[domain].size(); i++) {
		for (esint n = run::mesh->nodes->dintervals[domain][i].begin; n < run::mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = run::mesh->nodes->coordinates->datatarray()[n] - center;
			run::data->N1[domain].dense_values.push_back((   0 - r56 / r55 * ((-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y) - r46 / r44 * (-p.y / norm.x)) / norm.z);
			run::data->N1[domain].dense_values.push_back((-p.z - r56 / r55 * ((   0 - r45 / r44 * ( p.x / norm.x)) / norm.y) - r46 / r44 * ( p.x / norm.x)) / norm.z);
			run::data->N1[domain].dense_values.push_back(( p.y - r56 / r55 * (( p.x - r45 / r44 * (   0 / norm.x)) / norm.y) - r46 / r44 * (   0 / norm.x)) / norm.z);
		}
	}

	std::vector<esint> fixPoints;
//	if (_BEMDomain[domain]) {
//		fixPoints = std::vector<esint>(
//				run::mesh->FETIData->surfaceFixPoints.begin() + run::mesh->FETIData->sFixPointsDistribution[domain],
//				run::mesh->FETIData->surfaceFixPoints.begin() + run::mesh->FETIData->sFixPointsDistribution[domain + 1]);
//	} else {
		fixPoints = std::vector<esint>(
				run::mesh->FETIData->innerFixPoints.begin() + run::mesh->FETIData->iFixPointsDistribution[domain],
				run::mesh->FETIData->innerFixPoints.begin() + run::mesh->FETIData->iFixPointsDistribution[domain + 1]);
//	}

	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 6;
	Nt.cols = run::data->K[domain].cols;
	Nt.nnz  = 9 * fixPoints.size();
	Nt.type = 'G';

	std::vector<esint> &ROWS = Nt.CSR_I_row_indices;
	std::vector<esint> &COLS = Nt.CSR_J_col_indices;
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

	auto n2DOF = [&] (esint node) {
		auto dit = run::mesh->nodes->dintervals[domain].begin();
		while (dit->end < node) { ++dit; }
		return dit->DOFOffset + node - dit->begin;
	};

	for (size_t c = 0; c < 3; c++) {
		for (size_t i = 0; i < fixPoints.size(); i++) {
			COLS.push_back(3 * n2DOF(fixPoints[i]) + c + 1);
		}
	}
	VALS.insert(VALS.end(), 3 * fixPoints.size(), 1);

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = run::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 0 + 1);
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 1 + 1);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = run::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 0 + 1);
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 2 + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = run::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 1 + 1);
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 2 + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.y);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );

	run::data->RegMat[domain].MatMat(Nt, 'N', N);
	run::data->RegMat[domain].MatTranspose();
	run::data->RegMat[domain].RemoveLower();

	SparseSolverCPU NtN;
	NtN.ImportMatrix(run::data->RegMat[domain]);
	run::data->RegMat[domain].Clear();

	NtN.Factorization("Create RegMat");
	NtN.SolveMat_Sparse(Nt);
	NtN.Clear();

	run::data->RegMat[domain].MatMat(N, 'N', Nt);
	run::data->RegMat[domain].MatScale(run::data->K[domain].getDiagonalMaximum());

}




