
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HYPREPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HYPREPROVIDER_H_

#include "../provider.h"

namespace espreso {

class HYPREProvider: public Provider {

public:
	HYPREProvider(LoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;
	double& solutionPrecision();
};


}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HYPREPROVIDER_H_ */
