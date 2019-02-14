
#include "structuralmechanics.mklpdssprovider.h"

#include "config/ecf/physics/structuralmechanics.h"
#include "basis/matrices/matrixtype.h"

using namespace espreso;

StructuralMechanicsMKLPDSSProvider::StructuralMechanicsMKLPDSSProvider(DataHolder *data, StructuralMechanicsLoadStepConfiguration &configuration)
: MKLPDSSProvider(data, configuration), _configuration(configuration)
{

}

MatrixType StructuralMechanicsMKLPDSSProvider::getMatrixType() const
{
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}




