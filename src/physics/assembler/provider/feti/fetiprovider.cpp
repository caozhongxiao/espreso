
#include "fetiprovider.h"

#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "physics/assembler/dataholder.h"

#include "basis/logging/logging.h"
#include "basis/containers/serializededata.h"
#include "config/ecf/solver/feti.h"
#include "config/ecf/physics/physicssolver/loadstep.h"

#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/fetidatastore.h"

#include "solver/generic/SparseMatrix.h"


using namespace espreso;

FETIProvider::FETIProvider(DataHolder *data, LoadStepConfiguration &configuration)
: Provider(data, configuration)
{
	_data->N1.clear();
	_data->N2.clear();
	_data->RegMat.clear();

	_data->B0.clear();
	_data->B0subdomainsMap.clear();

	_data->N1.resize(info::mesh->elements->ndomains);
	_data->N2.resize(info::mesh->elements->ndomains);
	_data->RegMat.resize(info::mesh->elements->ndomains);

	_data->B0.resize(info::mesh->elements->ndomains);
	_data->B0subdomainsMap.resize(info::mesh->elements->ndomains);

	if (_configuration.type == LoadStepConfiguration::TYPE::TRANSIENT) {
		_data->computeKernelCallback = [&] (FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster) {};
		_data->computeKernelsCallback = [&] (FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster) {};

		_data->computeKernelsFromOrigKCallback = [&] (FETI_REGULARIZATION regularization, size_t scSize, bool ortogonalCluster) {
			_data->K.swap(_data->origK);
			_data->N1.swap(_data->origKN1);
			_data->N2.swap(_data->origKN2);
			_data->RegMat.swap(_data->origRegMat);
			makeStiffnessMatricesRegular(regularization, scSize, ortogonalCluster);

			_data->K.swap(_data->origK);
			_data->N1.swap(_data->origKN1);
			_data->N2.swap(_data->origKN2);
			_data->RegMat.swap(_data->origRegMat);
		};

		_data->computeKernelFromOrigKCallback = [&] (FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster) {
			_data->K[domain].swap(_data->origK[domain]);
			_data->N1[domain].swap(_data->origKN1[domain]);
			_data->N2[domain].swap(_data->origKN2[domain]);
			_data->RegMat[domain].swap(_data->origRegMat[domain]);
			makeStiffnessMatrixRegular(regularization, scSize, domain, ortogonalCluster);

			_data->K[domain].swap(_data->origK[domain]);
			_data->N1[domain].swap(_data->origKN1[domain]);
			_data->N2[domain].swap(_data->origKN2[domain]);
			_data->RegMat[domain].swap(_data->origRegMat[domain]);
		};
	} else {
		_data->computeKernelsCallback = [&] (FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster) {
			makeStiffnessMatricesRegular(regularization, scSize, ortogonalCluster);
		};

		_data->computeKernelCallback = [&] (FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster) {
			makeStiffnessMatrixRegular(regularization, scSize, domain, ortogonalCluster);
		};
	}

	_data->assembleB0Callback = [&] (FETI_B0_TYPE type, const std::vector<SparseMatrix> &kernels) {
		switch (type) {
		case FETI_B0_TYPE::CORNERS:
			assembleB0FromCorners();
			break;
		case FETI_B0_TYPE::KERNELS:
			assembleB0FromKernels(kernels);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown type of B0";
		}
	};
}

bool FETIProvider::needOriginalStiffnessMatrices()
{
	// TODO: some solvers may need in other cases
	return false || Provider::needOriginalStiffnessMatrices();
}

double& FETIProvider::solutionPrecision()
{
	return _configuration.feti.precision;
}
void FETIProvider::makeStiffnessMatricesRegular(FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster)
{
	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		makeStiffnessMatrixRegular(regularization, scSize, d, ortogonalCluster);
		ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	ESINFO(PROGRESS3);
}

void FETIProvider::makeStiffnessMatrixRegular(FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster)
{
	switch (regularization) {

	case FETI_REGULARIZATION::ANALYTIC:
		analyticRegularization(domain, ortogonalCluster);
		_data->RegMat[domain].RemoveLower();
		_data->K[domain].MatAddInPlace(_data->RegMat[domain], 'N', 1);
		_data->RegMat[domain].ConvertToCOO(1);
		break;

	case FETI_REGULARIZATION::ALGEBRAIC:
		switch (_data->K[domain].mtype) {
			double norm;
			esint defect;

		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			_data->K[domain].get_kernel_from_K(_data->K[domain], _data->RegMat[domain], _data->N1[domain], norm, defect, domain, scSize);
			break;

		case MatrixType::REAL_UNSYMMETRIC:
			_data->K[domain].get_kernels_from_nonsym_K(_data->K[domain], _data->RegMat[domain], _data->N1[domain], _data->N2[domain], norm, defect, domain, scSize);
			break;

		default:
			ESINFO(ERROR) << "Unknown matrix type for regularization.";
		}
		break;
	}
}

void FETIProvider::assembleUniformB0FromCorners(int DOFs)
{
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		_data->B0[d].cols = _data->K[d].cols;
	}

	size_t lambdas = 1;

	auto domains = info::mesh->FETIData->cornerDomains->cbegin();
	for (size_t n = 0; n < info::mesh->FETIData->corners.size(); ++n, ++domains) {
		for (size_t dof = 0; dof < 1; dof++) {
			for (size_t d1 = 0, d2 = 1; d2 < domains->size(); ++d1, ++d2) {

				auto d1it = info::mesh->nodes->dintervals[domains->at(d1)].begin();
				auto d2it = info::mesh->nodes->dintervals[domains->at(d2)].begin();
				while (d1it->end < info::mesh->FETIData->corners[n]) {
					++d1it;
				}
				while (d2it->end < info::mesh->FETIData->corners[n]) {
					++d2it;
				}

				for (int dof = 0; dof < DOFs; dof++) {
					_data->B0[domains->at(d1)].I_row_indices.push_back(lambdas);
					_data->B0[domains->at(d1)].J_col_indices.push_back(DOFs * (info::mesh->FETIData->corners[n] - d1it->begin + d1it->DOFOffset) + dof + 1);
					_data->B0[domains->at(d1)].V_values.push_back(1);

					_data->B0[domains->at(d2)].I_row_indices.push_back(lambdas);
					_data->B0[domains->at(d2)].J_col_indices.push_back(DOFs * (info::mesh->FETIData->corners[n] - d2it->begin + d2it->DOFOffset) + dof + 1);
					_data->B0[domains->at(d2)].V_values.push_back(-1);

					lambdas++;
				}
			}
		}
	}

	#pragma omp parallel for
	for  (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_data->B0[p].rows = lambdas - 1;
		_data->B0[p].cols = _data->K[p].cols;
		_data->B0[p].nnz = _data->B0[p].I_row_indices.size();

		_data->B0subdomainsMap[p].reserve(_data->B0[p].nnz);
		for (esint i = _data->B0subdomainsMap[p].size(); i < _data->B0[p].nnz; i++) {
			_data->B0subdomainsMap[p].push_back(_data->B0[p].I_row_indices[i] - 1);
		}
	}
}

void FETIProvider::assembleUniformB0FromKernels(const std::vector<SparseMatrix> &kernels, int DOFs)
{
	std::vector<esint> rowIndex(info::mesh->FETIData->inodesDomains.size());
	std::vector<esint> rCounters(*std::max_element(info::mesh->elements->clusters.begin(), info::mesh->elements->clusters.end()) + 1);

	for (size_t i = 0; i < info::mesh->FETIData->inodesDomains.size(); i++) {
		esint domain = info::mesh->FETIData->inodesDomains[i].first;
		esint ndomain = info::mesh->FETIData->inodesDomains[i].second;
		esint cluster = info::mesh->elements->clusters[domain];
		rowIndex[i] = rCounters[cluster];
		rCounters[cluster] += std::max(kernels[domain].cols, std::max(kernels[ndomain].cols, (esint)1));
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &nodes = info::mesh->FETIData->interfaceNodes->datatarray();
		int sign, cols, master;
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			for (size_t i = 0; i < info::mesh->FETIData->inodesDomains.size(); i++) {
				sign = 0;
				if (info::mesh->FETIData->inodesDomains[i].first == d) {
					sign = 1;
				}
				if (info::mesh->FETIData->inodesDomains[i].second == d) {
					sign = -1;
				}

				if (sign != 0) {
					master = info::mesh->FETIData->inodesDomains[i].first;
					if (kernels[master].cols < kernels[info::mesh->FETIData->inodesDomains[i].second].cols) {
						master = info::mesh->FETIData->inodesDomains[i].second;
					}
					cols = kernels[master].cols;
					if (cols) {
						for (esint c = 0; c < cols; c++) {
							auto dit = info::mesh->nodes->dintervals[d].begin();
							auto masterit = info::mesh->nodes->dintervals[master].begin();
							for (esint n = info::mesh->FETIData->inodesDistribution[i]; n < info::mesh->FETIData->inodesDistribution[i + 1]; n++) {
								while (dit->end < nodes[n]) {
									++dit;
								}
								while (masterit->end < nodes[n]) {
									++masterit;
								}
								for (int dof = 0; dof < DOFs; dof++) {
									_data->B0[d].I_row_indices.push_back(rowIndex[i] + c + 1);
									_data->B0[d].J_col_indices.push_back(DOFs * (dit->DOFOffset + nodes[n] - dit->begin) + dof + 1);
									_data->B0[d].V_values.push_back(sign * kernels[master].dense_values[kernels[master].rows * c + DOFs * (masterit->DOFOffset + nodes[n] - masterit->begin) + dof]);
								}
							}
						}
					} else {
						auto dit = info::mesh->nodes->dintervals[d].begin();
						for (esint n = info::mesh->FETIData->inodesDistribution[i]; n < info::mesh->FETIData->inodesDistribution[i + 1]; n++) {
							while (dit->end < nodes[n]) {
								++dit;
							}
							for (int dof = 0; dof < DOFs; dof++) {
								_data->B0[d].I_row_indices.push_back(rowIndex[i] + 1);
								_data->B0[d].J_col_indices.push_back(DOFs * (dit->DOFOffset + nodes[n] - dit->begin) + dof + 1);
								_data->B0[d].V_values.push_back(sign);
							}
						}
					}
				}
			}
			_data->B0[d].rows = rCounters[info::mesh->elements->clusters[d]];
			_data->B0[d].cols = _data->K[d].cols;
			_data->B0[d].nnz = _data->B0[d].I_row_indices.size();
			_data->B0subdomainsMap[d].reserve(_data->B0[d].nnz);
			for (esint i = _data->B0subdomainsMap[d].size(); i < _data->B0[d].nnz; i++) {
				_data->B0subdomainsMap[d].push_back(_data->B0[d].I_row_indices[i]);
			}
		}
	}

}




