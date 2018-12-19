
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_GLOBAL_HEATTRANSFER_GLOBALPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_GLOBAL_HEATTRANSFER_GLOBALPROVIDER_H_

#include "globalprovider.h"

namespace espreso {

struct HeatTransferLoadStepConfiguration;

class HeatTransferGlobalProvider: public GlobalProvider {

public:
	HeatTransferGlobalProvider(HeatTransferLoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;

protected:
	HeatTransferLoadStepConfiguration &_configuration;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_GLOBAL_HEATTRANSFER_GLOBALPROVIDER_H_ */
