
#include "heattransfer.mklpdssprovider.h"

#include "config/ecf/physics/heattransfer.h"
#include "basis/matrices/matrixtype.h"

using namespace espreso;

HeatTransferMKLPDSSProvider::HeatTransferMKLPDSSProvider(DataHolder *data, HeatTransferLoadStepConfiguration &configuration)
: MKLPDSSProvider(data, configuration), _configuration(configuration)
{

}

MatrixType HeatTransferMKLPDSSProvider::getMatrixType() const
{
	if (	(_configuration.translation_motions.size()) ||
			(_configuration.mode == LoadStepConfiguration::MODE::NONLINEAR && _configuration.nonlinear_solver.tangent_matrix_correction)) {

		return MatrixType::REAL_UNSYMMETRIC;
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}




