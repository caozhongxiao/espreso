
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HEATTRANSFER_HYPREPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HEATTRANSFER_HYPREPROVIDER_H_

#include "hypreprovider.h"

namespace espreso {

struct HeatTransferLoadStepConfiguration;

class HeatTransferHYPREProvider: public HYPREProvider {

public:
	HeatTransferHYPREProvider(HeatTransferLoadStepConfiguration &configuration);

protected:
	HeatTransferLoadStepConfiguration &_configuration;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_HYPRE_HEATTRANSFER_HYPREPROVIDER_H_ */
