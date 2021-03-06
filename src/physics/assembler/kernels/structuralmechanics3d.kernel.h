
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS3D_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS3D_KERNEL_H_

#include "structuralmechanics.kernel.h"

namespace espreso {

struct DenseMatrix;
struct Element;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
enum Matrices: int;
struct SolverParameters;

struct StructuralMechanics3DKernel: public StructuralMechanicsKernel
{
	struct ElementIterator {
		Element *element;

		double *coordinates;
		double *displacement;
		double *acceleration;
		double *angularVelocity;
		double *initialTemperature;
		double *temperature;

		const MaterialConfiguration *material;

		bool largeDisplacement;

		ElementIterator()
		: element(NULL),
		  coordinates(NULL),
		  displacement(NULL),
		  acceleration(NULL),
		  angularVelocity(NULL),
		  initialTemperature(NULL), temperature(NULL),
		  material(NULL),
		  largeDisplacement(NULL) {}
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

	void processElement(Matrices matrices, const SolverParameters &parameters, const ElementIterator &iterator, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Ce, DenseMatrix &Re, DenseMatrix &fe) const;
	void processFace(Matrices matrices, const SolverParameters &parameters, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const;
	void processEdge(Matrices matrices, const SolverParameters &parameters, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const;

	void processSolution(const SolutionIterator &iterator);

protected:
	void assembleLinearElasticMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, DenseMatrix &K) const;
	void assembleHyperElasticMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, DenseMatrix &K) const;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_STRUCTURALMECHANICS3D_KERNEL_H_ */
