
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_STRUCTURALMECHANICS3D_FETIPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_STRUCTURALMECHANICS3D_FETIPROVIDER_H_

#include "structuralmechanics.fetiprovider.h"

namespace espreso {

struct StructuralMechanicsLoadStepConfiguration;

class StructuralMechanics3DFETIProvider: public StructuralMechanicsFETIProvider {

public:
	StructuralMechanics3DFETIProvider(StructuralMechanicsLoadStepConfiguration &configuration);

protected:
	void analyticRegularization(esint domain, bool ortogonalCluster);

	void assembleB0FromCorners()
	{
		assembleUniformB0FromCorners(3);
	}

	void assembleB0FromKernels(const std::vector<SparseMatrix> &kernels)
	{
		assembleUniformB0FromKernels(kernels, 3);
	}
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_STRUCTURALMECHANICS3D_FETIPROVIDER_H_ */
