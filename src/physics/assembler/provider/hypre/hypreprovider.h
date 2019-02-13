
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HYPREPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HYPREPROVIDER_H_

#include "physics/assembler/provider/provider.h"

namespace espreso {

class HYPREProvider: public Provider {

public:
	HYPREProvider(DataHolder *data, LoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;

	virtual bool needMatrixVectorProduct();
	virtual bool needOriginalStiffnessMatrices();

	double& solutionPrecision();
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HYPREPROVIDER_H_ */
