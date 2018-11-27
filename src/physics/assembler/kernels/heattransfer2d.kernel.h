
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER2D_KERNEL_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER2D_KERNEL_H_

#include "heattransfer.kernel.h"

namespace espreso {

struct DenseMatrix;
struct Element;
struct Step;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
enum Matrices: int;

struct HeatTransfer2DKernel: public HeatTransferKernel
{
	struct ElementIterator {
		Element *element;

		double *temperature;
		double *coordinates;
		double *gradient;
		double *motion;
		double *heat;
		double *thickness;

		const MaterialConfiguration *material;

		ElementIterator()
		: element(NULL),
		  temperature(NULL), coordinates(NULL),
		  gradient(NULL), motion(NULL), heat(NULL), thickness(NULL),
		  material(NULL) {}
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

	HeatTransfer2DKernel(const HeatTransferGlobalSettings &settings, const HeatTransferOutputSettings &output);

	void processElement(Matrices matrices, const ElementIterator &iterator, const Step &step, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
	void processEdge(Matrices matrices, const BoundaryIterator &iterator, const Step &step, DenseMatrix &Ke, DenseMatrix &fe) const;
//	void processNode(eslocal domain, const BoundaryRegionStore *region, Matrices matrices, eslocal nindex, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe) const;
//	void processSolution();

protected:
	void assembleMaterialMatrix(eslocal node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const;
//	void postProcessElement(eslocal domain, eslocal eindex);
};

}



#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_HEATTRANSFER2D_KERNEL_H_ */
