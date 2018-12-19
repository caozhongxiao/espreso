
#include "fetiprovider.h"

#include "../../../../globals/run.h"
#include "../../../../basis/logging/logging.h"
#include "../../../../basis/containers/serializededata.h"
#include "../../../../config/ecf/solver/feti.h"
#include "../../../../config/ecf/physics/physicssolver/loadstep.h"

#include "../../../dataholder.h"

#include "../../../../mesh/mesh.h"
#include "../../../../mesh/store/elementstore.h"
#include "../../../../mesh/store/nodestore.h"
#include "../../../../mesh/store/fetidatastore.h"

#include "../../../../solver/generic/SparseMatrix.h"


using namespace espreso;

FETIProvider::FETIProvider(LoadStepConfiguration &configuration)
: _configuration(configuration)
{

}

bool FETIProvider::needOriginalStiffnessMatrices()
{
	return _configuration.type == LoadStepConfiguration::TYPE::TRANSIENT;
}

void FETIProvider::makeStiffnessMatricesRegular(FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster)
{
	#pragma omp parallel for
	for (eslocal d = 0; d < run::mesh->elements->ndomains; d++) {
		makeStiffnessMatrixRegular(regularization, scSize, d, ortogonalCluster);
		ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	ESINFO(PROGRESS3);
}

void FETIProvider::makeStiffnessMatrixRegular(FETI_REGULARIZATION regularization, int scSize, eslocal domain, bool ortogonalCluster)
{
	switch (regularization) {

	case FETI_REGULARIZATION::ANALYTIC:
		analyticRegularization(domain, ortogonalCluster);
		run::data->RegMat[domain].RemoveLower();
		run::data->K[domain].MatAddInPlace(run::data->RegMat[domain], 'N', 1);
		run::data->RegMat[domain].ConvertToCOO(1);
		break;

	case FETI_REGULARIZATION::ALGEBRAIC:
		switch (run::data->K[domain].mtype) {
			double norm;
			eslocal defect;

		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			run::data->K[domain].get_kernel_from_K(run::data->K[domain], run::data->RegMat[domain], run::data->N1[domain], norm, defect, domain, scSize);
			break;

		case MatrixType::REAL_UNSYMMETRIC:
			run::data->K[domain].get_kernels_from_nonsym_K(run::data->K[domain], run::data->RegMat[domain], run::data->N1[domain], run::data->N2[domain], norm, defect, domain, scSize);
			break;

		default:
			ESINFO(ERROR) << "Unknown matrix type for regularization.";
		}
		break;
	}
}

void FETIProvider::assembleB0FromCorners(int DOFs)
{
	for (size_t d = 0; d < run::mesh->elements->ndomains; d++) {
		run::data->B0[d].cols = run::data->K[d].cols;
	}

	size_t lambdas = 1;

	auto domains = run::mesh->FETIData->cornerDomains->cbegin();
	for (size_t n = 0; n < run::mesh->FETIData->corners.size(); ++n, ++domains) {
		for (size_t dof = 0; dof < 1; dof++) {
			for (size_t d1 = 0, d2 = 1; d2 < domains->size(); ++d1, ++d2) {

				auto d1it = run::mesh->nodes->dintervals[domains->at(d1)].begin();
				auto d2it = run::mesh->nodes->dintervals[domains->at(d2)].begin();
				while (d1it->end < run::mesh->FETIData->corners[n]) {
					++d1it;
				}
				while (d2it->end < run::mesh->FETIData->corners[n]) {
					++d2it;
				}

				for (int dof = 0; dof < DOFs; dof++) {
					run::data->B0[domains->at(d1)].I_row_indices.push_back(lambdas);
					run::data->B0[domains->at(d1)].J_col_indices.push_back(DOFs * (run::mesh->FETIData->corners[n] - d1it->begin + d1it->DOFOffset) + dof + 1);
					run::data->B0[domains->at(d1)].V_values.push_back(1);

					run::data->B0[domains->at(d2)].I_row_indices.push_back(lambdas);
					run::data->B0[domains->at(d2)].J_col_indices.push_back(DOFs * (run::mesh->FETIData->corners[n] - d2it->begin + d2it->DOFOffset) + dof + 1);
					run::data->B0[domains->at(d2)].V_values.push_back(-1);

					lambdas++;
				}
			}
		}
	}

	#pragma omp parallel for
	for  (size_t p = 0; p < run::mesh->elements->ndomains; p++) {
		run::data->B0[p].rows = lambdas - 1;
		run::data->B0[p].cols = run::data->K[p].cols;
		run::data->B0[p].nnz = run::data->B0[p].I_row_indices.size();

		run::data->B0subdomainsMap[p].reserve(run::data->B0[p].nnz);
		for (eslocal i = run::data->B0subdomainsMap[p].size(); i < run::data->B0[p].nnz; i++) {
			run::data->B0subdomainsMap[p].push_back(run::data->B0[p].I_row_indices[i] - 1);
		}
	}
}

void FETIProvider::assembleB0FromKernels(const std::vector<SparseMatrix> &kernels, int DOFs)
{
	std::vector<eslocal> rowIndex(run::mesh->FETIData->inodesDomains.size());
	std::vector<eslocal> rCounters(*std::max_element(run::mesh->elements->clusters.begin(), run::mesh->elements->clusters.end()) + 1);

	for (size_t i = 0; i < run::mesh->FETIData->inodesDomains.size(); i++) {
		eslocal domain = run::mesh->FETIData->inodesDomains[i].first;
		eslocal ndomain = run::mesh->FETIData->inodesDomains[i].second;
		eslocal cluster = run::mesh->elements->clusters[domain];
		rowIndex[i] = rCounters[cluster];
		rCounters[cluster] += std::max(kernels[domain].cols, std::max(kernels[ndomain].cols, (eslocal)1));
	}

	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &nodes = run::mesh->FETIData->interfaceNodes->datatarray();
		int sign, cols, master;
		for (eslocal d = run::mesh->elements->domainDistribution[t]; d < run::mesh->elements->domainDistribution[t + 1]; d++) {
			for (size_t i = 0; i < run::mesh->FETIData->inodesDomains.size(); i++) {
				sign = 0;
				if (run::mesh->FETIData->inodesDomains[i].first == d) {
					sign = 1;
				}
				if (run::mesh->FETIData->inodesDomains[i].second == d) {
					sign = -1;
				}

				if (sign != 0) {
					master = run::mesh->FETIData->inodesDomains[i].first;
					if (kernels[master].cols < kernels[run::mesh->FETIData->inodesDomains[i].second].cols) {
						master = run::mesh->FETIData->inodesDomains[i].second;
					}
					cols = kernels[master].cols;
					if (cols) {
						for (eslocal c = 0; c < cols; c++) {
							auto dit = run::mesh->nodes->dintervals[d].begin();
							auto masterit = run::mesh->nodes->dintervals[master].begin();
							for (eslocal n = run::mesh->FETIData->inodesDistribution[i]; n < run::mesh->FETIData->inodesDistribution[i + 1]; n++) {
								while (dit->end < nodes[n]) {
									++dit;
								}
								while (masterit->end < nodes[n]) {
									++masterit;
								}
								for (int dof = 0; dof < DOFs; dof++) {
									run::data->B0[d].I_row_indices.push_back(rowIndex[i] + c + 1);
									run::data->B0[d].J_col_indices.push_back(DOFs * (dit->DOFOffset + nodes[n] - dit->begin) + dof + 1);
									run::data->B0[d].V_values.push_back(sign * kernels[master].dense_values[kernels[master].rows * c + DOFs * (masterit->DOFOffset + nodes[n] - masterit->begin) + dof]);
								}
							}
						}
					} else {
						auto dit = run::mesh->nodes->dintervals[d].begin();
						for (eslocal n = run::mesh->FETIData->inodesDistribution[i]; n < run::mesh->FETIData->inodesDistribution[i + 1]; n++) {
							while (dit->end < nodes[n]) {
								++dit;
							}
							for (int dof = 0; dof < DOFs; dof++) {
								run::data->B0[d].I_row_indices.push_back(rowIndex[i] + 1);
								run::data->B0[d].J_col_indices.push_back(DOFs * (dit->DOFOffset + nodes[n] - dit->begin) + dof + 1);
								run::data->B0[d].V_values.push_back(sign);
							}
						}
					}
				}
			}
			run::data->B0[d].rows = rCounters[run::mesh->elements->clusters[d]];
			run::data->B0[d].cols = run::data->K[d].cols;
			run::data->B0[d].nnz = run::data->B0[d].I_row_indices.size();
			run::data->B0subdomainsMap[d].reserve(run::data->B0[d].nnz);
			for (eslocal i = run::data->B0subdomainsMap[d].size(); i < run::data->B0[d].nnz; i++) {
				run::data->B0subdomainsMap[d].push_back(run::data->B0[d].I_row_indices[i]);
			}
		}
	}

}




