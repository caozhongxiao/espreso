
#include "heattransfer.hypreprovider.h"

#include "../../../../config/ecf/physics/heattransfer.h"

using namespace espreso;

HeatTransferHYPREProvider::HeatTransferHYPREProvider(HeatTransferLoadStepConfiguration &configuration)
: HYPREProvider(configuration), _configuration(configuration)
{

}
