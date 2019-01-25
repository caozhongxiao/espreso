
#include "hypreprovider.h"

#include "basis/matrices/matrixtype.h"
#include "config/ecf/physics/structuralmechanics.h"

using namespace espreso;

HYPREProvider::HYPREProvider(DataHolder *data, LoadStepConfiguration &configuration)
: Provider(data, configuration)
{

}

MatrixType HYPREProvider::getMatrixType() const
{
	return MatrixType::REAL_UNSYMMETRIC; // Hypre always get full matrix
}

bool HYPREProvider::needMatrixVectorProduct()
{
	return Provider::needMatrixVectorProduct() || _configuration.mode == LoadStepConfiguration::MODE::NONLINEAR;
}

bool HYPREProvider::needOriginalStiffnessMatrices()
{
	return true;
//	return Provider::needOriginalStiffnessMatrices() || _configuration.mode == LoadStepConfiguration::MODE::NONLINEAR;
}

double& HYPREProvider::solutionPrecision()
{
	return _configuration.multigrid.precision;
}





