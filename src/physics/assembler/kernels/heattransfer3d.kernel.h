
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER3D_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER3D_KERNEL_H_

#include "heattransfer.kernel.h"

namespace espreso {

struct DenseMatrix;
struct Element;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
enum Matrices: int;

struct HeatTransfer3DKernel: public HeatTransferKernel
{
	struct ElementIterator {
		Element *element;

		double *temperature;
		double *coordinates;
		double *gradient;
		double *motion;
		double *heat;

		const MaterialConfiguration *material;

		ElementIterator()
		: element(NULL),
		  temperature(NULL), coordinates(NULL),
		  gradient(NULL), motion(NULL), heat(NULL),
		  material(NULL) {}
	};

	struct SolutionIterator: public ElementIterator {
		double *phase;
		double *latentHeat;
		double *gradient;
		double *flux;

		SolutionIterator()
		: ElementIterator(),
		  phase(NULL), latentHeat(NULL),
		  gradient(NULL), flux(NULL) {}
	};

	struct BoundaryIterator {
		Element *element;

		double *temperature;
		double *coordinates;
		double *thickness;

		double *htc;
		double *emissivity;
		double *externalTemperature;
		double *heatFlow;
		double *heatFlux;

		double regionArea;

		bool radiation, convection;

		BoundaryIterator()
		: element(NULL),
		  temperature(NULL), coordinates(NULL), thickness(NULL),
		  htc(NULL), emissivity(NULL), externalTemperature(NULL),
		  heatFlow(NULL), heatFlux(NULL),
		  regionArea(0),
		  radiation(false), convection(false) {}

	};

	HeatTransfer3DKernel(const HeatTransferGlobalSettings &settings, const HeatTransferOutputSettings &output);

	void processElement(Matrices matrices, const ElementIterator &iterator, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processFace(Matrices matrices, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const;
	void processEdge(Matrices matrices, const BoundaryIterator &iterator, DenseMatrix &Ke, DenseMatrix &fe) const;

	void processSolution(const SolutionIterator &iterator);

protected:
	void assembleMaterialMatrix(eslocal node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const;

};

}


#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER3D_KERNEL_H_ */
