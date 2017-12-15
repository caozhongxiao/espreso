
#ifndef SRC_ASSEMBLER_PHYSICS_LAPLACESTEKLOVPOINCARE3D_H_
#define SRC_ASSEMBLER_PHYSICS_LAPLACESTEKLOVPOINCARE3D_H_

#include "heattransfer3d.h"

namespace espreso {

struct LaplaceSteklovPoincare3D: public HeatTransfer3D
{
	LaplaceSteklovPoincare3D(Mesh *mesh, Instance *instance, Step *step, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	void prepareHybridTotalFETIWithKernels();

	virtual void preprocessData();

	virtual void updateMatrix(Matrices matrices, size_t domain);
	virtual void updateMatrix(Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe);

	void processSolution();
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_LAPLACESTEKLOVPOINCARE3D_H_ */
