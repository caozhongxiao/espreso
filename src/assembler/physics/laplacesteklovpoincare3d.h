
#ifndef SRC_ASSEMBLER_PHYSICS_LAPLACESTEKLOVPOINCARE3D_H_
#define SRC_ASSEMBLER_PHYSICS_LAPLACESTEKLOVPOINCARE3D_H_

#include "heattransfer3d.h"

namespace espreso {

struct LaplaceSteklovPoincare3D: public HeatTransfer3D
{
	LaplaceSteklovPoincare3D(Mesh *mesh, Instance *instance, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	void prepareHybridTotalFETIWithKernels();

	virtual void preprocessData(const Step &step);

	virtual void updateMatrix(const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution);
	virtual void updateMatrix(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution);

	void processSolution(const Step &step);

protected:
	enum BEMSolutionIndex: size_t {
		BOUNDARY = 0,

		SIZE     = 1
	};

	static size_t BEMOffset;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_LAPLACESTEKLOVPOINCARE3D_H_ */
