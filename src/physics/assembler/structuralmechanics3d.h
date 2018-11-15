
#ifndef SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS3D_H_
#define SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS3D_H_

#include "../assembler/physics3d.h"
#include "../assembler/structuralmechanics.h"

namespace espreso {

class MaterialBaseConfiguration;

struct StructuralMechanics3D: public StructuralMechanics, public Physics3D
{
	StructuralMechanics3D(Mesh *mesh, Instance *instance, Step *step, const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	virtual void initLocalDOFs(std::vector<eslocal> &offsets) { initLocalNodeUniformDOFs(offsets, 3); }
	virtual void initGlobalDOFs(eslocal &offset) { initGlobalNodeUniformDOFs(offset, 3); }

	virtual void buildLocalCSRPattern() { buildLocalNodeUniformCSRPattern(3); }
	virtual void buildGlobalCSRPattern() { buildGlobalNodeUniformCSRPattern(3); }

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


#endif /* SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS3D_H_ */
