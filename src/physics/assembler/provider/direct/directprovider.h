
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_DIRECT_DIRECTPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_DIRECT_DIRECTPROVIDER_H_

#include "physics/assembler/provider/provider.h"

namespace espreso {

class DirectProvider: public Provider {

public:
	DirectProvider(DataHolder *data, LoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;

	double& solutionPrecision() { return _precision; }

protected:
	double _precision;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_DIRECT_DIRECTPROVIDER_H_ */
