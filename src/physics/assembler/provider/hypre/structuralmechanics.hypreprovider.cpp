
#include "structuralmechanics.hypreprovider.h"

#include "config/ecf/physics/structuralmechanics.h"

using namespace espreso;

StructuralMechanicsHYPREProvider::StructuralMechanicsHYPREProvider(StructuralMechanicsLoadStepConfiguration &configuration)
: HYPREProvider(configuration), _configuration(configuration)
{

}
