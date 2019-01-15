
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_STRUCTURALMECHANICS2D_FETIPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_STRUCTURALMECHANICS2D_FETIPROVIDER_H_

#include "structuralmechanics.fetiprovider.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanics2DFETIProvider: public StructuralMechanicsFETIProvider {

public:
	StructuralMechanics2DFETIProvider(DataHolder *data, StructuralMechanicsLoadStepConfiguration &configuration);

protected:
	void analyticRegularization(esint domain, bool ortogonalCluster);

	void assembleB0FromCorners()
	{
		assembleUniformB0FromCorners(2);
	}

	void assembleB0FromKernels(const std::vector<SparseMatrix> &kernels)
	{
		assembleUniformB0FromKernels(kernels, 2);
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_STRUCTURALMECHANICS2D_FETIPROVIDER_H_ */
