
#include "provider.h"

#include "../../../config/ecf/physics/physicssolver/loadstep.h"

using namespace espreso;

Provider::Provider(LoadStepConfiguration &configuration)
: _configuration(configuration)
{

}

bool Provider::needOriginalStiffnessMatrices()
{
	return _configuration.type == LoadStepConfiguration::TYPE::TRANSIENT;
}



