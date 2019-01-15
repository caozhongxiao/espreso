
#include "heattransfer.fetiprovider.h"

#include "physics/assembler/dataholder.h"
#include "physics/assembler/composer/feti/feticomposer.h"

#include "esinfo/meshinfo.h"

#include "basis/logging/logging.h"
#include "basis/matrices/matrixtype.h"
#include "config/ecf/physics/heattransfer.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "solver/generic/SparseMatrix.h"


using namespace espreso;


HeatTransferFETIProvider::HeatTransferFETIProvider(DataHolder *data, HeatTransferLoadStepConfiguration &configuration)
: FETIProvider(data, configuration), _configuration(configuration)
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
			ElementsRegionStore *region = info::mesh->eregion(it->first);
			for (esint i = info::mesh->elements->eintervalsDistribution[domain]; i < info::mesh->elements->eintervalsDistribution[domain + 1]; i++) {
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
	if (_data->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC";
	}

	if (
			(_configuration.feti.conjugate_projector != FETI_CONJ_PROJECTOR::CONJ_R &&
			_configuration.feti.conjugate_projector != FETI_CONJ_PROJECTOR::CONJ_K) ||
			_configuration.type != LoadStepConfiguration::TYPE::TRANSIENT) {

		if (_configuration.convection.size()) {
			for (auto it = _configuration.convection.begin(); it != _configuration.convection.end(); ++it) {
				BoundaryRegionStore *region = info::mesh->bregion(it->first);
				if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
					return;
				}
			}
		}

		if (_configuration.diffuse_radiation.size()) {
			for (auto it = _configuration.diffuse_radiation.begin(); it != _configuration.diffuse_radiation.end(); ++it) {
				BoundaryRegionStore *region = info::mesh->bregion(it->first);
				if (region->eintervalsDistribution[domain] != region->eintervalsDistribution[domain + 1]) {
					return;
				}
			}
		}
	}

	if (_configuration.feti.conjugate_projector == FETI_CONJ_PROJECTOR::CONJ_K) {
		_data->N1[domain].rows = _data->K[domain].rows;
		_data->N1[domain].cols = 1;
		_data->N1[domain].nnz = _data->N1[domain].rows * _data->N1[domain].cols;
		_data->N1[domain].type = 'G';

		std::vector<double> diagonal = _data->K[domain].getDiagonal();
		_data->N1[domain].dense_values.insert(_data->N1[domain].dense_values.end(), diagonal.begin(), diagonal.end());

		_data->RegMat[domain].rows = _data->K[domain].rows;
		_data->RegMat[domain].cols = _data->K[domain].cols;
		_data->RegMat[domain].nnz  = 1;
		_data->RegMat[domain].type = _data->K[domain].type;

		_data->RegMat[domain].I_row_indices.push_back(1);
		_data->RegMat[domain].J_col_indices.push_back(1);
		_data->RegMat[domain].V_values.push_back(_data->K[domain].getDiagonalMaximum());
		_data->RegMat[domain].ConvertToCSR(1);
	} else {
		double value;
		if (ortogonalCluster) {
			size_t nSum = 0;
			for (size_t d = 0; d < info::mesh->elements->ndomains; d++) {
				if (info::mesh->elements->clusters[d] == info::mesh->elements->clusters[domain]) {
					nSum += _data->K[d].rows;
				}
			}
			value = 1 / sqrt(nSum);
		} else {
			value = 1 / sqrt(_data->K[domain].rows);
		}

		_data->N1[domain].rows = _data->K[domain].rows;
		_data->N1[domain].cols = 1;
		_data->N1[domain].nnz = _data->N1[domain].rows * _data->N1[domain].cols;
		_data->N1[domain].type = 'G';

		_data->N1[domain].dense_values.resize(_data->N1[domain].nnz, value);

		_data->RegMat[domain].rows = _data->K[domain].rows;
		_data->RegMat[domain].cols = _data->K[domain].cols;
		_data->RegMat[domain].nnz  = 1;
		_data->RegMat[domain].type = _data->K[domain].type;

		_data->RegMat[domain].I_row_indices.push_back(1);
		_data->RegMat[domain].J_col_indices.push_back(1);
		_data->RegMat[domain].V_values.push_back(_data->K[domain].getDiagonalMaximum());
		_data->RegMat[domain].ConvertToCSR(1);
	}
}

