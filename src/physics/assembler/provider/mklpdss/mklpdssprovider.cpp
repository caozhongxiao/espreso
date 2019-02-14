
#include "mklpdssprovider.h"

using namespace espreso;

MKLPDSSProvider::MKLPDSSProvider(DataHolder *data, LoadStepConfiguration &configuration)
: Provider(data, configuration), _precision(1e-5)
{

}


