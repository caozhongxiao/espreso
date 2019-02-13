
#include "directprovider.h"
#include "basis/matrices/matrixtype.h"

using namespace espreso;

DirectProvider::DirectProvider(DataHolder *data, LoadStepConfiguration &configuration)
: Provider(data, configuration), _precision(1e-5)
{

}

MatrixType DirectProvider::getMatrixType() const
{
	return MatrixType::REAL_UNSYMMETRIC;
}


