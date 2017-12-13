
#ifndef SRC_ASSEMBLER_PHYSICS_HEATTRANSFER3D_H_
#define SRC_ASSEMBLER_PHYSICS_HEATTRANSFER3D_H_

#include "heattransfer.h"
#include "physics3d.h"

namespace espreso {

class MaterialBaseConfiguration;

struct HeatTransfer3D: public HeatTransfer, public Physics3D
{
	HeatTransfer3D(Mesh *mesh, Instance *instance, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	virtual void processElement(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	virtual void processFace(const Step &step, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	virtual void processEdge(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	virtual void processNode(const Step &step, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	virtual void processSolution(const Step &step);

protected:
	void assembleMaterialMatrix(const Step &step, eslocal eindex, eslocal node, const Point &p, const MaterialBaseConfiguration *mat, double phase, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const;
	void postProcessElement(const Step &step, eslocal eindex);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_HEATTRANSFER3D_H_ */
