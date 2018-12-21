
#include "structuralmechanics.hypreprovider.h"

#include "../../../../basis/matrices/matrixtype.h"
#include "../../../../config/ecf/physics/structuralmechanics.h"

using namespace espreso;

StructuralMechanicsHYPREProvider::StructuralMechanicsHYPREProvider(StructuralMechanicsLoadStepConfiguration &configuration)
: Provider(configuration), _configuration(configuration)
{

}

MatrixType StructuralMechanicsHYPREProvider::getMatrixType() const
{
	return MatrixType::REAL_UNSYMMETRIC; // Hypre always get full matrix
}
