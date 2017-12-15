
#ifndef SRC_ASSEMBLER_PHYSICS_LAMESTEKLOVPOINCARE3D_H_
#define SRC_ASSEMBLER_PHYSICS_LAMESTEKLOVPOINCARE3D_H_

#include "structuralmechanics3d.h"

namespace espreso {

struct LameSteklovPoincare3D: public StructuralMechanics3D
{
	LameSteklovPoincare3D(Mesh *mesh, Instance *instance, Step *step, const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	void prepareHybridTotalFETIWithKernels();

	virtual void preprocessData();

	virtual void updateMatrix(Matrices matrices, size_t domain);
	virtual void updateMatrix(Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe);

	void processSolution();
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LAMESTEKLOVPOINCARE3D_H_ */
