
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_GLOBAL_STRUCTURALMECHANICS_GLOBALPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_GLOBAL_STRUCTURALMECHANICS_GLOBALPROVIDER_H_

#include "globalprovider.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsGlobalProvider: public GlobalProvider {

public:
	StructuralMechanicsGlobalProvider(StructuralMechanicsLoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;

protected:
	StructuralMechanicsLoadStepConfiguration &_configuration;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_GLOBAL_STRUCTURALMECHANICS_GLOBALPROVIDER_H_ */
