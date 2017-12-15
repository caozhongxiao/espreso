
#ifndef SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS3D_H_
#define SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS3D_H_

#include "physics3d.h"
#include "structuralmechanics.h"

namespace espreso {

struct StructuralMechanics3D: public StructuralMechanics, public Physics3D
{
	StructuralMechanics3D(Mesh *mesh, Instance *instance, Step *step, const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	void analyticRegularization(size_t domain, bool ortogonalCluster);

	void processElement(eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processSolution();

protected:
	void assembleMaterialMatrix(eslocal eindex, eslocal node, double temp, DenseMatrix &K) const;
	void postProcessElement(eslocal eindex);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS3D_H_ */
