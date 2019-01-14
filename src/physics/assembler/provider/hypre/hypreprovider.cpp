
#include "hypreprovider.h"

#include "basis/matrices/matrixtype.h"
#include "config/ecf/physics/structuralmechanics.h"

using namespace espreso;

HYPREProvider::HYPREProvider(LoadStepConfiguration &configuration)
: Provider(configuration)
{

}

MatrixType HYPREProvider::getMatrixType() const
{
	return MatrixType::REAL_UNSYMMETRIC; // Hypre always get full matrix
}

double& HYPREProvider::solutionPrecision()
{
	return _configuration.multigrid.precision;
}





