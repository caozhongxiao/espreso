
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS3D_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS3D_KERNEL_H_

#include "structuralmechanics.kernel.h"

namespace espreso {

struct DenseMatrix;
struct Element;
struct Step;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
enum Matrices: int;

struct StructuralMechanics3DKernel: public StructuralMechanicsKernel
{
	struct ElementIterator {
		Element *element;

		double *coordinates;
		double *acceleration;
		double *initialTemperature;
		double *temperature;

		const MaterialConfiguration *material;

		ElementIterator()
		: element(NULL),
		  coordinates(NULL),
		  acceleration(NULL),
		  initialTemperature(NULL), temperature(NULL),
		  material(NULL) {}
	};

	struct SolutionIterator: public ElementIterator {
		SolutionIterator()
		: ElementIterator() {}
	};

	struct BoundaryIterator {
		Element *element;

		double *coordinates;
		double *normalPressure;

		BoundaryIterator()
		: element(NULL),
		  coordinates(NULL),
		  normalPressure(NULL) {}

	};

	StructuralMechanics3DKernel(const StructuralMechanicsGlobalSettings &settings, const StructuralMechanicsOutputSettings &output);

	void processElement(Matrices matrices, const ElementIterator &iterator, const Step &step, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processFace(Matrices matrices, const BoundaryIterator &iterator, const Step &step, DenseMatrix &Ke, DenseMatrix &fe) const;
	void processEdge(Matrices matrices, const BoundaryIterator &iterator, const Step &step, DenseMatrix &Ke, DenseMatrix &fe) const;

	void processSolution(const SolutionIterator &iterator, const Step &step);

protected:
	void assembleMaterialMatrix(eslocal node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, DenseMatrix &K) const;

};

}


#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS3D_KERNEL_H_ */
