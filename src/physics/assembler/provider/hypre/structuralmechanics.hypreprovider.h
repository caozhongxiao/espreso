
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_STRUCTURALMECHANICS_HYPREPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_STRUCTURALMECHANICS_HYPREPROVIDER_H_

#include "../provider.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsHYPREProvider: public Provider {

public:
	StructuralMechanicsHYPREProvider(StructuralMechanicsLoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;

protected:
	StructuralMechanicsLoadStepConfiguration &_configuration;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_STRUCTURALMECHANICS_HYPREPROVIDER_H_ */
