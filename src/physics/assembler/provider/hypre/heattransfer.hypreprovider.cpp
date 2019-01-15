
#include "heattransfer.hypreprovider.h"

#include "config/ecf/physics/heattransfer.h"

using namespace espreso;

HeatTransferHYPREProvider::HeatTransferHYPREProvider(DataHolder *data, HeatTransferLoadStepConfiguration &configuration)
: HYPREProvider(data, configuration), _configuration(configuration)
{

}
