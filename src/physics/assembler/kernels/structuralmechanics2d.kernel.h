
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS2D_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS2D_KERNEL_H_

#include "structuralmechanics.kernel.h"

namespace espreso {

struct DenseMatrix;
struct Element;
struct Step;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
enum Matrices: int;

struct StructuralMechanics2DKernel: public StructuralMechanicsKernel
{
	struct ElementIterator {
		Element *element;

		double *coordinates;
		double *acceleration;
		double *angularVelocity;
		double *initialTemperature;
		double *temperature;
		double *thickness;

		const MaterialConfiguration *material;

		ElementIterator()
		: element(NULL),
		  coordinates(NULL),
		  acceleration(NULL), angularVelocity(NULL),
		  initialTemperature(NULL), temperature(NULL),
		  thickness(NULL),
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
		double *thickness;

		BoundaryIterator()
		: element(NULL),
		  coordinates(NULL),
		  normalPressure(NULL), thickness(NULL) {}

	};

	StructuralMechanics2DKernel(const StructuralMechanicsGlobalSettings &settings, const StructuralMechanicsOutputSettings &output);

	void processElement(Matrices matrices, const ElementIterator &iterator, const Step &step, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processEdge(Matrices matrices, const BoundaryIterator &iterator, const Step &step, DenseMatrix &Ke, DenseMatrix &fe) const;

	void processSolution(const SolutionIterator &iterator, const Step &step);

protected:
	void assembleMaterialMatrix(eslocal node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, DenseMatrix &K) const;

};

}



#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS2D_KERNEL_H_ */
