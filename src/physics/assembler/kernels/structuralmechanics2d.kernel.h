
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS2D_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS2D_KERNEL_H_

#include "structuralmechanics.kernel.h"

namespace espreso {

struct DenseMatrix;
struct Element;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
enum Matrices: int;
struct SolverParameters;

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

	void processElement(Matrices matrices, const SolverParameters &parameters, const ElementIterator &iterator, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processEdge(Matrices matrices, const SolverParameters &parameters, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const;

	void processSolution(const SolutionIterator &iterator);

protected:
	void assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, DenseMatrix &K) const;

};

}



#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS2D_KERNEL_H_ */
