
#ifndef SRC_ASSEMBLER_PHYSICS_HEATTRANSFER2D_H_
#define SRC_ASSEMBLER_PHYSICS_HEATTRANSFER2D_H_

#include "heattransfer.h"
#include "physics2d.h"

namespace espreso {

enum class Property;
class MaterialBaseConfiguration;

struct HeatTransfer2D: public HeatTransfer, public Physics2D
{
	HeatTransfer2D(Mesh *mesh, Instance *instance, Step *step, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	void processElement(eslocal domain, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processFace(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processEdge(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processSolution();

protected:
	void assembleMaterialMatrix(eslocal eindex, eslocal node, const Point &p, const MaterialBaseConfiguration *mat, double phase, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const;
	void postProcessElement(eslocal domain, eslocal eindex);
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_HEATTRANSFER2D_H_ */
