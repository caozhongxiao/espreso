
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_MKLPDSS_STRUCTURALMECHANICS_MKLPDSSPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_MKLPDSS_STRUCTURALMECHANICS_MKLPDSSPROVIDER_H_

#include "mklpdssprovider.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanicsMKLPDSSProvider: public MKLPDSSProvider {

public:
	StructuralMechanicsMKLPDSSProvider(DataHolder *data, StructuralMechanicsLoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;

protected:
	StructuralMechanicsLoadStepConfiguration &_configuration;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_MKLPDSS_STRUCTURALMECHANICS_MKLPDSSPROVIDER_H_ */
