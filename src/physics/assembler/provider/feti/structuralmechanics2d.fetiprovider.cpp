
#include "structuralmechanics2d.fetiprovider.h"

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


StructuralMechanics2DFETIProvider::StructuralMechanics2DFETIProvider(StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanicsFETIProvider(configuration)
{

}

void StructuralMechanics2DFETIProvider::analyticRegularization(eslocal domain, bool ortogonalCluster)
{
	Point center; size_t np; double norm;
	if (ortogonalCluster) {
		center = _cCenter[run::mesh->elements->clusters[domain]];
		np = _cNp[run::mesh->elements->clusters[domain]];
		norm = _cNorm[run::mesh->elements->clusters[domain]].x;
	} else {
		center = _dCenter[domain];
		np = _dNp[domain];
		norm = _dNorm[domain].x;
	}

	run::data->N1[domain].rows = run::data->K[domain].rows;
	run::data->N1[domain].cols = 3;
	run::data->N1[domain].nnz = run::data->N1[domain].rows * run::data->N1[domain].cols;
	run::data->N1[domain].type = 'G';

	run::data->N1[domain].dense_values.reserve(run::data->N1[domain].nnz);

	for (size_t c = 0; c < 2; c++) {
		std::vector<double> kernel = { 0, 0 };
		kernel[c] = 1 / std::sqrt(np);
		for (size_t i = 0; i < run::data->K[domain].rows / 2; i++) {
			run::data->N1[domain].dense_values.insert(run::data->N1[domain].dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < run::mesh->nodes->dintervals[domain].size(); i++) {
		for (eslocal n = run::mesh->nodes->dintervals[domain][i].begin; n < run::mesh->nodes->dintervals[domain][i].end; ++n) {
			Point p = run::mesh->nodes->coordinates->datatarray()[n] - center;
			run::data->N1[domain].dense_values.push_back(-p.y / norm);
			run::data->N1[domain].dense_values.push_back( p.x / norm);
		}
	}

	std::vector<eslocal> fixPoints;
//	if (_BEMDomain[domain]) {
//		fixPoints = std::vector<eslocal>(
//				run::mesh->FETIData->surfaceFixPoints.begin() + run::mesh->FETIData->sFixPointsDistribution[domain],
//				run::mesh->FETIData->surfaceFixPoints.begin() + run::mesh->FETIData->sFixPointsDistribution[domain + 1]);
//	} else {
		fixPoints = std::vector<eslocal>(
				run::mesh->FETIData->innerFixPoints.begin() + run::mesh->FETIData->iFixPointsDistribution[domain],
				run::mesh->FETIData->innerFixPoints.begin() + run::mesh->FETIData->iFixPointsDistribution[domain + 1]);
//	}

	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 3;
	Nt.cols = run::data->K[domain].cols;
	Nt.nnz  = 4 * fixPoints.size();
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
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());

	auto n2DOF = [&] (eslocal node) {
		auto dit = run::mesh->nodes->dintervals[domain].begin();
		while (dit->end < node) { ++dit; }
		return dit->DOFOffset + node - dit->begin;
	};

	for (size_t c = 0; c < 2; c++) {
		for (size_t i = 0; i < fixPoints.size(); i++) {
			COLS.push_back(2 * n2DOF(fixPoints[i]) + c + 1);
		}
	}
	VALS.insert(VALS.end(), 2 * fixPoints.size(), 1);

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = run::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
		COLS.push_back(2 * n2DOF(fixPoints[i]) + 0 + 1);
		COLS.push_back(2 * n2DOF(fixPoints[i]) + 1 + 1);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
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




