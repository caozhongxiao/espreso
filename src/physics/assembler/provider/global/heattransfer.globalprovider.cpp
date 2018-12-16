
#include "heattransfer.globalprovider.h"

#include "../../../../basis/matrices/matrixtype.h"
#include "../../../../config/ecf/physics/heattransfer.h"

using namespace espreso;

HeatTransferGlobalProvider::HeatTransferGlobalProvider(HeatTransferLoadStepConfiguration &configuration)
: _configuration(configuration)
{

}

MatrixType HeatTransferGlobalProvider::getMatrixType() const
{
	if (_configuration.translation_motions.size()) {
		return MatrixType::REAL_UNSYMMETRIC;
	}
	return MatrixType::REAL_UNSYMMETRIC; // Hypre always get full matrix
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}
