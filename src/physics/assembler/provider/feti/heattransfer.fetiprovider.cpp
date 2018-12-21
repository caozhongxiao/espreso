
#include "heattransfer.fetiprovider.h"

#include "../../../dataholder.h"
#include "../../../../globals/run.h"
#include "../../../../basis/logging/logging.h"
#include "../../../../basis/matrices/matrixtype.h"
#include "../../../../config/ecf/physics/heattransfer.h"

#include "../../../../mesh/mesh.h"
#include "../../../../mesh/store/elementstore.h"
#include "../../../../mesh/store/elementsregionstore.h"
#include "../../../../solver/generic/SparseMatrix.h"


using namespace espreso;


HeatTransferFETIProvider::HeatTransferFETIProvider(HeatTransferLoadStepConfiguration &configuration)
: FETIProvider(configuration), _configuration(configuration)
{

}

MatrixType HeatTransferFETIProvider::getMatrixType() const
{
	if (	(_configuration.translation_motions.size()) ||
			(_configuration.mode == LoadStepConfiguration::MODE::NONLINEAR && _configuration.nonlinear_solver.tangent_matrix_correction)) {

		return MatrixType::REAL_UNSYMMETRIC;
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

MatrixType HeatTransferFETIProvider::getMatrixType(esint domain) const
{
	if (_configuration.mode == LoadStepConfiguration::MODE::NONLINEAR && _configuration.nonlinear_solver.tangent_matrix_correction) {
		return MatrixType::REAL_UNSYMMETRIC;
	}

	if (_configuration.translation_motions.size()) {
		for (auto it = _configuration.translation_motions.begin(); it != _configuration.translation_motions.end(); ++it) {
			ElementsRegionStore *region = run::mesh->eregion(it->first);
			for (esint i = run::mesh->elements->eintervalsDistribution[domain]; i < run::mesh->elements->eintervalsDistribution[domain + 1]; i++) {
				if (region->eintervals[i].begin != region->eintervals[i].end) {
					return MatrixType::REAL_UNSYMMETRIC;
				}
			}
		}
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

void HeatTransferFETIProvider::analyticRegularization(esint domain, bool ortogonalCluster)
{
	if (run::data->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC";
	}

	if (
			(_configuration.feti.conjugate_projector != FETI_CONJ_PROJECTOR::CONJ_R &&
			_configuration.feti.conjugate_projector != FETI_CONJ_PROJECTOR::CONJ_K) ||
			_configuration.type != LoadStepConfiguration::TYPE::TRANSIENT) {

		if (_configuration.convection.size()) {
			for (auto it = _configuration.convection.begin(); it != _configuration.convection.end(); ++it) {
				BoundaryRegionStore *region = run::mesh->bregion(it->first);
				if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
					return;
				}
			}
		}

		if (_configuration.diffuse_radiation.size()) {
			for (auto it = _configuration.diffuse_radiation.begin(); it != _configuration.diffuse_radiation.end(); ++it) {
				BoundaryRegionStore *region = run::mesh->bregion(it->first);
				if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
					return;
				}
			}
		}
	}

	if (_configuration.feti.conjugate_projector == FETI_CONJ_PROJECTOR::CONJ_K) {
		run::data->N1[domain].rows = run::data->K[domain].rows;
		run::data->N1[domain].cols = 1;
		run::data->N1[domain].nnz = run::data->N1[domain].rows * run::data->N1[domain].cols;
		run::data->N1[domain].type = 'G';

		std::vector<double> diagonal = run::data->K[domain].getDiagonal();
		run::data->N1[domain].dense_values.insert(run::data->N1[domain].dense_values.end(), diagonal.begin(), diagonal.end());

		run::data->RegMat[domain].rows = run::data->K[domain].rows;
		run::data->RegMat[domain].cols = run::data->K[domain].cols;
		run::data->RegMat[domain].nnz  = 1;
		run::data->RegMat[domain].type = run::data->K[domain].type;

		run::data->RegMat[domain].I_row_indices.push_back(1);
		run::data->RegMat[domain].J_col_indices.push_back(1);
		run::data->RegMat[domain].V_values.push_back(run::data->K[domain].getDiagonalMaximum());
		run::data->RegMat[domain].ConvertToCSR(1);
	} else {
		double value;
		if (ortogonalCluster) {
			size_t nSum = 0;
			for (size_t d = 0; d < run::mesh->elements->ndomains; d++) {
				if (run::mesh->elements->clusters[d] == run::mesh->elements->clusters[domain]) {
					nSum += run::data->K[d].rows;
				}
			}
			value = 1 / sqrt(nSum);
		} else {
			value = 1 / sqrt(run::data->K[domain].rows);
		}

		run::data->N1[domain].rows = run::data->K[domain].rows;
		run::data->N1[domain].cols = 1;
		run::data->N1[domain].nnz = run::data->N1[domain].rows * run::data->N1[domain].cols;
		run::data->N1[domain].type = 'G';

		run::data->N1[domain].dense_values.resize(run::data->N1[domain].nnz, value);

		run::data->RegMat[domain].rows = run::data->K[domain].rows;
		run::data->RegMat[domain].cols = run::data->K[domain].cols;
		run::data->RegMat[domain].nnz  = 1;
		run::data->RegMat[domain].type = run::data->K[domain].type;

		run::data->RegMat[domain].I_row_indices.push_back(1);
		run::data->RegMat[domain].J_col_indices.push_back(1);
		run::data->RegMat[domain].V_values.push_back(run::data->K[domain].getDiagonalMaximum());
		run::data->RegMat[domain].ConvertToCSR(1);
	}
}

