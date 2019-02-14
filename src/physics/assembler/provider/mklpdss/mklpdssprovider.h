
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_MKLPDSS_MKLPDSSPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_MKLPDSS_MKLPDSSPROVIDER_H_

#include "physics/assembler/provider/provider.h"

namespace espreso {

class MKLPDSSProvider: public Provider {

public:
	MKLPDSSProvider(DataHolder *data, LoadStepConfiguration &configuration);

	double& solutionPrecision() { return _precision; }

protected:
	double _precision; // dummy since MKL PDSS does not have this option
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_MKLPDSS_MKLPDSSPROVIDER_H_ */
