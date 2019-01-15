
#include "provider.h"

#include "config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

Provider::Provider(DataHolder *data, LoadStepConfiguration &configuration)
: _data(data), _configuration(configuration)
{

}

bool Provider::needOriginalStiffnessMatrices()
{
	return _configuration.type == LoadStepConfiguration::TYPE::TRANSIENT;
}



