
#ifndef SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_FETIPROVIDER_H_
#define SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_FETIPROVIDER_H_

#include "physics/assembler/provider/provider.h"
#include <vector>

namespace espreso {

enum class MatrixType;
enum class FETI_REGULARIZATION;
class SparseMatrix;

class FETIProvider: public Provider {

public:
	FETIProvider(LoadStepConfiguration &configuration);

	virtual MatrixType getMatrixType() const =0;
	virtual MatrixType getMatrixType(esint domain) const =0;

	bool needOriginalStiffnessMatrices();
	double& solutionPrecision();

protected:
	void makeStiffnessMatricesRegular(FETI_REGULARIZATION regularization, int scSize, bool ortogonalCluster);
	void makeStiffnessMatrixRegular(FETI_REGULARIZATION regularization, int scSize, esint domain, bool ortogonalCluster);

	virtual void assembleB0FromCorners() =0;
	virtual void assembleB0FromKernels(const std::vector<SparseMatrix> &kernels) =0;

	void assembleUniformB0FromCorners(int DOFs);
	void assembleUniformB0FromKernels(const std::vector<SparseMatrix> &kernels, int DOFs);

	virtual void analyticRegularization(esint domain, bool ortogonalCluster) =0;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_PROVIDER_FETI_FETIPROVIDER_H_ */
