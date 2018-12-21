
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_STRUCTURALMECHANICS_HYPREPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_STRUCTURALMECHANICS_HYPREPROVIDER_H_

#include "hypreprovider.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsHYPREProvider: public HYPREProvider {

public:
	StructuralMechanicsHYPREProvider(StructuralMechanicsLoadStepConfiguration &configuration);

protected:
	StructuralMechanicsLoadStepConfiguration &_configuration;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_STRUCTURALMECHANICS_HYPREPROVIDER_H_ */
