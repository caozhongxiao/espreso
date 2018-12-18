
#include "fetiprovider.h"

#include "../../../../globals/run.h"
#include "../../../../basis/logging/logging.h"
#include "../../../../config/ecf/solver/feti.h"
#include "../../../../config/ecf/physics/physicssolver/loadstep.h"

#include "../../../dataholder.h"

#include "../../../../mesh/mesh.h"
#include "../../../../mesh/store/elementstore.h"

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




