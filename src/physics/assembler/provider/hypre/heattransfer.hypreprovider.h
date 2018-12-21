
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HEATTRANSFER_HYPREPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HEATTRANSFER_HYPREPROVIDER_H_

#include "../provider.h"

namespace espreso {

struct HeatTransferLoadStepConfiguration;

class HeatTransferHYPREProvider: public Provider {

public:
	HeatTransferHYPREProvider(HeatTransferLoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;

protected:
	HeatTransferLoadStepConfiguration &_configuration;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HEATTRANSFER_HYPREPROVIDER_H_ */
