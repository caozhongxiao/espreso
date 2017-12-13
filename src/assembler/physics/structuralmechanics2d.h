
#ifndef SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS2D_H_
#define SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS2D_H_

#include "physics2d.h"
#include "structuralmechanics.h"

namespace espreso {

struct StructuralMechanics2D: public StructuralMechanics, public Physics2D
{
	StructuralMechanics2D(Mesh *mesh, Instance *instance, const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	void analyticRegularization(size_t domain, bool ortogonalCluster);

	void processElement(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processFace(const Step &step, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processEdge(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processNode(const Step &step, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processSolution(const Step &step);

protected:
	void assembleMaterialMatrix(const Step &step, eslocal eindex, eslocal node, double temp, DenseMatrix &K) const;
	void postProcessElement(const Step &step, eslocal eindex);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS2D_H_ */
