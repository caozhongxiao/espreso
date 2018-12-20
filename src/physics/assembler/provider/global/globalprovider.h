
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_GLOBAL_GLOBALPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_GLOBAL_GLOBALPROVIDER_H_

#include "../provider.h"

namespace espreso {

enum class MatrixType;

class GlobalProvider: public Provider {

public:
	GlobalProvider(LoadStepConfiguration &configuration);

	virtual MatrixType getMatrixType() const =0;
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_GLOBAL_GLOBALPROVIDER_H_ */
