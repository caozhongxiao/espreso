
#ifndef SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS3D_H_
#define SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS3D_H_

#include "physics3d.h"
#include "structuralmechanics.h"

namespace espreso {

struct StructuralMechanics3D: public StructuralMechanics, public Physics3D
{
	StructuralMechanics3D(Mesh *mesh, Instance *instance, const StructuralMechanicsConfiguration &configuration, const ResultsSelectionConfiguration &propertiesConfiguration);

	virtual std::vector<std::pair<ElementType, Property> > propertiesToStore() const;

	void prepare();
	void analyticRegularization(size_t domain, bool ortogonalCluster);

	void processElement(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processFace(const Step &step, Matrices matrices, eslocal findex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processEdge(const Step &step, Matrices matrices, eslocal eindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processNode(const Step &step, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processSolution(const Step &step);

	const std::vector<Property>& pointDOFs() const
	{
		static std::vector<Property> pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };
		return pointDOFs;
	}
	const std::vector<Property>& midPointDOFs() const
	{
		static std::vector<Property> midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };
		return midPointDOFs;
	}

protected:
	void assembleMaterialMatrix(const Step &step, eslocal eindex, eslocal node, double temp, DenseMatrix &K) const;
	void postProcessElement(const Step &step, eslocal eindex, std::vector<Solution*> &solution);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS3D_H_ */
