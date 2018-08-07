
#ifndef SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS2D_H_
#define SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS2D_H_

#include "physics2d.h"
#include "structuralmechanics.h"

namespace espreso {

struct MaterialBaseConfiguration;

struct StructuralMechanics2D: public StructuralMechanics, public Physics2D
{
	StructuralMechanics2D(Mesh *mesh, Instance *instance, Step *step, const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	void analyticRegularization(size_t domain, bool ortogonalCluster);

	void processBEM(eslocal domain, Matrices matrices);
	void processElement(eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processSolution();

protected:
	void assembleMaterialMatrix(eslocal node, const Point &p, const MaterialBaseConfiguration *mat, double temp, DenseMatrix &K) const;
	void postProcessElement(eslocal eindex);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS2D_H_ */
