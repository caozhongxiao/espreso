
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
: _configuration(configuration)
{
	run::data->N1.clear();
	run::data->N2.clear();
	run::data->RegMat.clear();

	run::data->N1.resize(run::mesh->elements->ndomains);
	run::data->N2.resize(run::mesh->elements->ndomains);
	run::data->RegMat.resize(run::mesh->elements->ndomains);

	if (_configuration.type == LoadStepConfiguration::TYPE::TRANSIENT) {
		run::data->computeKernelCallback = [&] (FETI_REGULARIZATION regularization, int scSize, eslocal domain, bool ortogonalCluster) {};
		run::data->computeKernelsCallback = [&] (FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster) {};

		run::data->computeKernelsFromOrigKCallback = [&] (FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster) {
			run::data->K.swap(run::data->origK);
			run::data->N1.swap(run::data->origKN1);
			run::data->N2.swap(run::data->origKN2);
			run::data->RegMat.swap(run::data->origRegMat);
			makeStiffnessMatricesRegular(regularization, scSize, ortogonalCluster);

			run::data->K.swap(run::data->origK);
			run::data->N1.swap(run::data->origKN1);
			run::data->N2.swap(run::data->origKN2);
			run::data->RegMat.swap(run::data->origRegMat);
		};

		run::data->computeKernelFromOrigKCallback = [&] (FETI_REGULARIZATION regularization, int scSize, eslocal domain, bool ortogonalCluster) {
			run::data->K[domain].swap(run::data->origK[domain]);
			run::data->N1[domain].swap(run::data->origKN1[domain]);
			run::data->N2[domain].swap(run::data->origKN2[domain]);
			run::data->RegMat[domain].swap(run::data->origRegMat[domain]);
			makeStiffnessMatrixRegular(regularization, scSize, domain, ortogonalCluster);

			run::data->K[domain].swap(run::data->origK[domain]);
			run::data->N1[domain].swap(run::data->origKN1[domain]);
			run::data->N2[domain].swap(run::data->origKN2[domain]);
			run::data->RegMat[domain].swap(run::data->origRegMat[domain]);
		};
	} else {
		run::data->computeKernelsCallback = [&] (FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster) {
			makeStiffnessMatricesRegular(regularization, scSize, ortogonalCluster);
		};

		run::data->computeKernelCallback = [&] (FETI_REGULARIZATION regularization, int scSize, eslocal domain, bool ortogonalCluster) {
			makeStiffnessMatrixRegular(regularization, scSize, domain, ortogonalCluster);
		};
	}

	run::data->assembleB0Callback = [&] (FETI_B0_TYPE type, const std::vector<SparseMatrix> &kernels) {
		run::data->B0.clear();
		run::data->B0.resize(run::mesh->elements->ndomains);
		for (eslocal d = 0; d < run::mesh->elements->ndomains; d++) {
			run::data->B0[d].type = 'G';
			run::data->B0subdomainsMap[d].clear();
		}
		switch (type) {
		case FETI_B0_TYPE::CORNERS:
//				physics.assembleB0FromCorners();
			break;
		case FETI_B0_TYPE::KERNELS:
//				physics.assembleB0FromKernels(kernels);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown type of B0";
		}
	};
}

MatrixType HeatTransferFETIProvider::getMatrixType() const
{
	if (_configuration.translation_motions.size()) {
		return MatrixType::REAL_UNSYMMETRIC;
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

MatrixType HeatTransferFETIProvider::getMatrixType(eslocal domain) const
{
	//	if (_step.tangentMatrixCorrection) {
	//		return MatrixType::REAL_UNSYMMETRIC;
	//	}

	if (_configuration.translation_motions.size()) {
		for (auto it = _configuration.translation_motions.begin(); it != _configuration.translation_motions.end(); ++it) {
			ElementsRegionStore *region = run::mesh->eregion(it->first);
			for (eslocal i = run::mesh->elements->eintervalsDistribution[domain]; i < run::mesh->elements->eintervalsDistribution[domain + 1]; i++) {
				if (region->eintervals[i].begin != region->eintervals[i].end) {
					return MatrixType::REAL_UNSYMMETRIC;
				}
			}
		}
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

void HeatTransferFETIProvider::analyticRegularization(eslocal domain, bool ortogonalCluster)
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

