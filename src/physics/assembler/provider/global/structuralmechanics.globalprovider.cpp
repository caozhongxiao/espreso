
#include "structuralmechanics.globalprovider.h"

#include "../../../../basis/matrices/matrixtype.h"
#include "../../../../config/ecf/physics/structuralmechanics.h"

using namespace espreso;

StructuralMechanicsGlobalProvider::StructuralMechanicsGlobalProvider(StructuralMechanicsLoadStepConfiguration &configuration)
: _configuration(configuration)
{

}

MatrixType StructuralMechanicsGlobalProvider::getMatrixType() const
{
	return MatrixType::REAL_UNSYMMETRIC; // Hypre always get full matrix
}
