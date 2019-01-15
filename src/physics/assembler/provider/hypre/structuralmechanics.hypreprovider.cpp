
#include "structuralmechanics.hypreprovider.h"

#include "config/ecf/physics/structuralmechanics.h"

using namespace espreso;

StructuralMechanicsHYPREProvider::StructuralMechanicsHYPREProvider(DataHolder *data, StructuralMechanicsLoadStepConfiguration &configuration)
: HYPREProvider(data, configuration), _configuration(configuration)
{

}
