
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_HEATTRANSFER_FETIPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_HEATTRANSFER_FETIPROVIDER_H_

#include "fetiprovider.h"

namespace espreso {

struct HeatTransferLoadStepConfiguration;

class HeatTransferFETIProvider: public FETIProvider {

public:
	HeatTransferFETIProvider(HeatTransferLoadStepConfiguration &configuration);

	MatrixType getMatrixType() const;
	MatrixType getMatrixType(eslocal domain) const;

protected:
	void analyticRegularization(eslocal domain, bool ortogonalCluster);

	void assembleB0FromCorners()
	{
		assembleUniformB0FromCorners(1);
	}

	void assembleB0FromKernels(const std::vector<SparseMatrix> &kernels)
	{
		assembleUniformB0FromKernels(kernels, 1);
	}

	HeatTransferLoadStepConfiguration &_configuration;
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_HEATTRANSFER_FETIPROVIDER_H_ */
