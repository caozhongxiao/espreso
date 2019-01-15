
#include "esinfo/meshinfo.h"
#include "physics/assembler/dataholder.h"
#include "physics/assembler/composer/composer.h"
#include "structuralmechanics3d.fetiprovider.h"

#include "basis/containers/serializededata.h"
#include "basis/matrices/matrixtype.h"
#include "basis/logging/logging.h"
#include "config/ecf/physics/structuralmechanics.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/fetidatastore.h"

#include "solver/generic/SparseMatrix.h"
#include "solver/specific/sparsesolvers.h"

using namespace espreso;

StructuralMechanics3DFETIProvider::StructuralMechanics3DFETIProvider(DataHolder *data, StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanicsFETIProvider(data, configuration)
{

}

void StructuralMechanics3DFETIProvider::analyticRegularization(esint domain, bool ortogonalCluster)
{
	if (_data->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC";
	}

	Point center = _dCenter[domain], norm = _dNorm[domain];
	double r44 = _dr44[domain], r45 = _dr45[domain], r46 = _dr46[domain], r55 = _dr55[domain], r56 = _dr56[domain];
	size_t np = _dNp[domain];

	if (ortogonalCluster) {
		size_t cluster = info::mesh->elements->clusters[domain];
		center = _cCenter[cluster], norm = _cNorm[cluster];
		r44 = _cr44[cluster], r45 = _cr45[cluster], r46 = _cr46[cluster], r55 = _cr55[cluster], r56 = _cr56[cluster];
		np = _cNp[cluster];
	} else {
		center = _dCenter[domain], norm = _dNorm[domain];
		r44 = _dr44[domain], r45 = _dr45[domain], r46 = _dr46[domain], r55 = _dr55[domain], r56 = _dr56[domain];
		np = _dNp[domain];
	}

	_data->N1[domain].rows = _data->K[domain].rows;
	_data->N1[domain].cols = 6;
	_data->N1[domain].nnz = _data->N1[domain].rows * _data->N1[domain].cols;
	_data->N1[domain].type = 'G';

	_data->N1[domain].dense_values.reserve(_data->N1[domain].nnz);

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };
		kernel[c] = 1 / std::sqrt(np);
		for (size_t i = 0; i < _data->K[domain].rows / 3; i++) {
			_data->N1[domain].dense_values.insert(_data->N1[domain].dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < info::mesh->nodes->dintervals[domain].size(); i++) {
		for (esint n = info::mesh->nodes->dintervals[domain][i].begin; n < info::mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = info::mesh->nodes->coordinates->datatarray()[n] - center;
			_data->N1[domain].dense_values.push_back(-p.y / norm.x);
			_data->N1[domain].dense_values.push_back( p.x / norm.x);
			_data->N1[domain].dense_values.push_back(             0);
		}
	}

	for (size_t i = 0; i < info::mesh->nodes->dintervals[domain].size(); i++) {
		for (esint n = info::mesh->nodes->dintervals[domain][i].begin; n < info::mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = info::mesh->nodes->coordinates->datatarray()[n] - center;
			_data->N1[domain].dense_values.push_back((-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y);
			_data->N1[domain].dense_values.push_back((   0 - r45 / r44 * ( p.x / norm.x)) / norm.y);
			_data->N1[domain].dense_values.push_back(( p.x - r45 / r44 * (   0 / norm.x)) / norm.y);
		}
	}

	for (size_t i = 0; i < info::mesh->nodes->dintervals[domain].size(); i++) {
		for (esint n = info::mesh->nodes->dintervals[domain][i].begin; n < info::mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = info::mesh->nodes->coordinates->datatarray()[n] - center;
			_data->N1[domain].dense_values.push_back((   0 - r56 / r55 * ((-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y) - r46 / r44 * (-p.y / norm.x)) / norm.z);
			_data->N1[domain].dense_values.push_back((-p.z - r56 / r55 * ((   0 - r45 / r44 * ( p.x / norm.x)) / norm.y) - r46 / r44 * ( p.x / norm.x)) / norm.z);
			_data->N1[domain].dense_values.push_back(( p.y - r56 / r55 * (( p.x - r45 / r44 * (   0 / norm.x)) / norm.y) - r46 / r44 * (   0 / norm.x)) / norm.z);
		}
	}

	std::vector<esint> fixPoints;
//	if (_BEMDomain[domain]) {
//		fixPoints = std::vector<esint>(
//				info::mesh->FETIData->surfaceFixPoints.begin() + info::mesh->FETIData->sFixPointsDistribution[domain],
//				info::mesh->FETIData->surfaceFixPoints.begin() + info::mesh->FETIData->sFixPointsDistribution[domain + 1]);
//	} else {
		fixPoints = std::vector<esint>(
				info::mesh->FETIData->innerFixPoints.begin() + info::mesh->FETIData->iFixPointsDistribution[domain],
				info::mesh->FETIData->innerFixPoints.begin() + info::mesh->FETIData->iFixPointsDistribution[domain + 1]);
//	}

	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 6;
	Nt.cols = _data->K[domain].cols;
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
		auto dit = info::mesh->nodes->dintervals[domain].begin();
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
		const Point &p = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 0 + 1);
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 1 + 1);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 0 + 1);
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 2 + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 1 + 1);
		COLS.push_back(3 * n2DOF(fixPoints[i]) + 2 + 1);
		VALS.push_back(-p.z);
		VALS.push_back( p.y);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );

	_data->RegMat[domain].MatMat(Nt, 'N', N);
	_data->RegMat[domain].MatTranspose();
	_data->RegMat[domain].RemoveLower();

	SparseSolverCPU NtN;
	NtN.ImportMatrix(_data->RegMat[domain]);
	_data->RegMat[domain].Clear();

	NtN.Factorization("Create RegMat");
	NtN.SolveMat_Sparse(Nt);
	NtN.Clear();

	_data->RegMat[domain].MatMat(N, 'N', Nt);
	_data->RegMat[domain].MatScale(_data->K[domain].getDiagonalMaximum());

}




