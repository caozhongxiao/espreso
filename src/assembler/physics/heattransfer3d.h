
#ifndef SRC_ASSEMBLER_PHYSICS_HEATTRANSFER3D_H_
#define SRC_ASSEMBLER_PHYSICS_HEATTRANSFER3D_H_

#include "heattransfer.h"
#include "physics3d.h"

namespace espreso {

class MaterialBaseConfiguration;

struct HeatTransfer3D: public HeatTransfer, public Physics3D
{
	HeatTransfer3D(OldMesh *mesh, Instance *instance, const HeatTransferConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	virtual void prepare();

	virtual std::vector<std::pair<ElementType, Property> > propertiesToStore() const;

	virtual void processElement(const Step &step, Matrices matrices, const OldElement *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	virtual void processFace(const Step &step, Matrices matrices, const OldElement *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	virtual void processEdge(const Step &step, Matrices matrices, const OldElement *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	virtual void processNode(const Step &step, Matrices matrices, const OldElement *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	virtual void processSolution(const Step &step);

protected:
	void assembleMaterialMatrix(const Step &step, const OldElement *e, eslocal node, const MaterialBaseConfiguration *mat, double phase, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const;
	void postProcessElement(const Step &step, const OldElement *e, std::vector<Solution*> &solution);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_HEATTRANSFER3D_H_ */
