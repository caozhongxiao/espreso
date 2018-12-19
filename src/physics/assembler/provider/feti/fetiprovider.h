
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_FETIPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_FETIPROVIDER_H_

#include "../provider.h"
#include <vector>

namespace espreso {

enum class MatrixType;
enum class FETI_REGULARIZATION;
struct LoadStepConfiguration;
class SparseMatrix;

class FETIProvider: public Provider {

public:
	FETIProvider(LoadStepConfiguration &configuration);

	virtual MatrixType getMatrixType() const =0;
	virtual MatrixType getMatrixType(eslocal domain) const =0;

	bool needOriginalStiffnessMatrices();

protected:
	void makeStiffnessMatricesRegular(FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster);
	void makeStiffnessMatrixRegular(FETI_REGULARIZATION regularization, int scSize, eslocal domain, bool ortogonalCluster);

	virtual void assembleB0FromCorners() =0;
	virtual void assembleB0FromKernels(const std::vector<SparseMatrix> &kernels) =0;

	void assembleUniformB0FromCorners(int DOFs);
	void assembleUniformB0FromKernels(const std::vector<SparseMatrix> &kernels, int DOFs);

	virtual void analyticRegularization(eslocal domain, bool ortogonalCluster) =0;

	LoadStepConfiguration &_configuration;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_FETIPROVIDER_H_ */
