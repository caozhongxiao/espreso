
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_MKLPDSS_HEATTRANSFER_MKLPDSSPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_MKLPDSS_HEATTRANSFER_MKLPDSSPROVIDER_H_

#include "mklpdssprovider.h"

namespace espreso {

struct HeatTransferLoadStepConfiguration;

class HeatTransferMKLPDSSProvider: public MKLPDSSProvider {

public:
	HeatTransferMKLPDSSProvider(DataHolder *data, HeatTransferLoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;

protected:
	HeatTransferLoadStepConfiguration &_configuration;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_MKLPDSS_HEATTRANSFER_MKLPDSSPROVIDER_H_ */
