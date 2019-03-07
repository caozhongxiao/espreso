
#include "mklpdssprovider.h"

#include "config/ecf/physics/structuralmechanics.h"

using namespace espreso;

MKLPDSSProvider::MKLPDSSProvider(DataHolder *data, LoadStepSolverConfiguration &configuration)
: Provider(data, configuration), _precision(1e-5)
{

}

bool MKLPDSSProvider::needMatrixVectorProduct()
{
	return Provider::needMatrixVectorProduct() || _configuration.mode == LoadStepSolverConfiguration::MODE::NONLINEAR;
}

bool MKLPDSSProvider::needOriginalStiffnessMatrices()
{
	return true;
}
